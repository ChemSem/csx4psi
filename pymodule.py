#
#@BEGIN LICENSE
#
# csx4psi by Psi4 Developer, a plugin to:
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

import psi4
import re
import os
import inputparser
import math
import warnings
import driver
from molutil import *
import p4util
from p4util.exceptions import *


def writeCSX(name, **kwargs):
    """function to write the CSX file

    """

    if not psi4.get_global_option('WRITE_CSX'):
        return

    # import csx_api for csx writing
    import os
    import math
    import inspect
    #import openbabel
    import qcdb
    import qcdb.periodictable
    import csx2_api as api
    lowername = name.lower()
    # Make sure the molecule the user provided is the active one
    if ('molecule' in kwargs):
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = psi4.get_active_molecule()
    molecule.update_geometry()
    # Determine the derivative type
    calledby = inspect.stack()[1][3]
    derdict = {
        'energy': 0,
        'property': 0,
        'gradient': 1,
        'optimize': 1,
        'frequency': 2,
        'frequencies': 2,
        'hessian': 2,
        }
    dertype = derdict[calledby]
    hasFreq = False
    # Start to write the CSX file
    # First grab molecular information and energies from psi4
    geom = molecule.save_string_xyz()  # OB
    atomLine = geom.split('\n')  # OB

    # general molecular information
    atomNum = molecule.natom()
    molSym = molecule.schoenflies_symbol()
    molCharge = molecule.molecular_charge()
    molMulti = molecule.multiplicity()
    # energy information
    molBasis = psi4.get_global_option('BASIS')
    molSpin = psi4.get_global_option('REFERENCE')
    molMethod = psi4.get_global_option('WFN')
    mol1E = psi4.get_variable('ONE-ELECTRON ENERGY')
    mol2E = psi4.get_variable('TWO-ELECTRON ENERGY')
    molNE = psi4.get_variable('NUCLEAR REPULSION ENERGY')
    molPE = mol1E + mol2E
    molEE = psi4.get_variable('CURRENT ENERGY')
    # wavefunction information
    try:
        wfn = kwargs['wfn']
    except AttributeError:
        pass
    if wfn:
        molOrbE = wfn.epsilon_a()
        molOrbEb = wfn.epsilon_b()
        orbNmopi = wfn.nmopi()
        orbNsopi = wfn.nsopi()
        orbNum = wfn.nmo() if molOrbE else 0
        orbSNum = wfn.nso()
        molOrb = wfn.Ca()
        orbNirrep = wfn.nirrep()
        orbAotoso = wfn.aotoso()
        orbDoccpi = wfn.doccpi()
        orbSoccpi = wfn.soccpi()
        basisNbf = wfn.basisset().nbf()
        basisDim = psi4.Dimension(1, 'basisDim')
        basisDim.__setitem__(0, basisNbf)
        wfnRestricted = True
        orbE = []
        hlist = []
        orblist = []
        orbOcc = []
        molOrbmo = psi4.Matrix('molOrbmo', basisDim, orbNmopi)
        molOrbmo.gemm(False, False, 1.0, orbAotoso, molOrb, 0.0)
        if molSpin == 'UHF':
            wfnRestricted = False
            orbEb = []
            hlistCb = []
            orblistCb = []
            orbOccCb = []
            molOrbCb = wfn.Cb()
            molOrbmoCb = psi4.Matrix('molOrbmoCb', basisDim, orbNmopi)
            molOrbmoCb.gemm(False, False, 1.0, orbAotoso, molOrbCb, 0.0)
        count = 0
        eleExtra = 1 if wfnRestricted else 0
        for ih in range(orbNirrep):
            for iorb in range(orbNmopi.__getitem__(ih)):
                hlist.append(ih)
                orblist.append(iorb)
                if molOrbE:
                    orbE.append(molOrbE.get(count))
                eleNum = 1 if iorb < (orbDoccpi.__getitem__(ih) + orbSoccpi.__getitem__(ih)) else 0
                eleNum += eleExtra if iorb < orbDoccpi.__getitem__(ih) else 0
                orbOcc.append(eleNum)
                count += 1
        orbMos = sorted(zip(orbE, zip(hlist, orblist)))
        orbOccString = ' '.join(str(x) for x in sorted(orbOcc, reverse=True))
        orbCaString = []
        for imos in range(orbNum):
            (h, s) = orbMos[imos][1]
            orbCa = []
            for iso in range(orbSNum):
                orbEle = molOrbmo.get(h, iso, s)
                orbCa.append(orbEle)
            orbCaString.append(' '.join(str(x) for x in orbCa))
        orbEString = ' '.join(str(x) for x in sorted(orbE))
        # now for beta spin
        if not wfnRestricted:
            count = 0
            for ih in range(orbNirrep):
                for iorb in range(orbNmopi.__getitem__(ih)):
                    hlistCb.append(ih)
                    orblist.append(iorb)
                    if molOrbEb:
                        orbEb.append(molOrbEb.get(count))
                    eleNum = 1 if iorb < (orbDoccpi.__getitem__(ih) + orbSoccpi.__getitem__(ih)) else 0
                    if iorb < orbDoccpi.__getitem__(ih):
                        eleNum += eleExtra
                    orbOccCb.append(eleNum)
                    count += 1
            orbMosCb = sorted(zip(orbEb, zip(hlist, orblist)))
            orbOccCbString = ' '.join(str(x) for x in sorted(orbOccCb, reverse=True))
            orbCbString = []
            for imos in range(orbNum):
                (h, s) = orbMosCb[imos][1]
                orbCb = []
                for iso in range(orbSNum):
                    orbEle = molOrbmoCb.get(h, iso, s)
                    orbCb.append(orbEle)
                orbCbString.append(' '.join(str(x) for x in orbCb))
            orbEbString = ' '.join(str(x) for x in sorted(orbEb))
        #   orbColString = ' '.join(str(x) for x in orbCol)
        if wfnRestricted:
            wfn1 = api.waveFunctionType(
                orbitalCount=orbNum,
                orbitalOccupancies=orbOccString)
            orbe1 = api.stringArrayType(unit='gc:hartree')
            orbe1.set_valueOf_(orbEString)
            orbs1 = api.orbitalsType()
            for iorb in range(orbNum):
                orb1 = api.stringArrayType(id=iorb+1)
                orb1.set_valueOf_(orbCaString[iorb])
                orbs1.add_orbital(orb1)
            wfn1.set_orbitals(orbs1)
            wfn1.set_orbitalEnergies(orbe1)
        else:
            wfn1 = api.waveFunctionType(orbitalCount=orbNum)
            # alpha electron: 1.5
            orbe1 = api.stringArrayType(unit='gc:hartree')
            orbe1.set_valueOf_(orbEString)
            wfn1.set_alphaOrbitalEnergies(orbe1)
            wfn1.set_alphaOrbitalOccupancies(orbOccString)
            aorbs1 = api.orbitalsType()
            for iorb in range(orbNum):
                orb1 = api.stringArrayType(id=iorb+1)
                orb1.set_valueOf_(orbCaString[iorb])
                aorbs1.add_orbital(orb1)
            wfn1.set_alphaOrbitals(aorbs1)
            # beta electron: 1.5
            orbeb1 = api.stringArrayType(unit='gc:hartree')
            orbeb1.set_valueOf_(orbEbString)
            wfn1.set_betaOrbitalEnergies(orbeb1)
            wfn1.set_betaOrbitalOccupancies(orbOccCbString)
            borbs1 = api.orbitalsType()
            for iorb in range(orbNum):
                orb1 = api.stringArrayType(id=iorb+1)
                orb1.set_valueOf_(orbCbString[iorb])
                borbs1.add_orbital(orb1)
            wfn1.set_betaOrbitals(borbs1)
    # frequency information
    if dertype == 2:
        hasFreq = True
        molFreq = psi4.get_frequencies()
        molFreqNum = molFreq.dim(0)
        frq = []
        irInt = []
        for ifrq in range(molFreqNum):
            frq.append(molFreq.get(ifrq))
            irInt.append(0.0)
        frqString = ' '.join(str(x) for x in frq)
        intString = ' '.join(str(x) for x in irInt)
        normMod = psi4.get_normalmodes()
        normMdString = []
        count = 0
        for ifrq in range(molFreqNum):
            normM = []
            for iatm in range(atomNum):
                for ixyz in range(3):
                    normM.append(normMod.get(count))
                    count += 1
            normMdString.append(' '.join(str(x) for x in normM))
        vib1 = api.vibAnalysisType(vibrationCount=molFreqNum)
        freq1 = api.stringArrayType(unit="gc:cm-1")
        freq1.set_valueOf_(frqString)
        vib1.set_frequencies(freq1)
        irint1 = api.stringArrayType()
        irint1.set_valueOf_(intString)
        vib1.set_irIntensities(irint1)
        norms1 = api.normalModesType()
        for ifrq in range(molFreqNum):
            norm1 = api.normalModeType(id=ifrq+1)
            norm1.set_valueOf_(normMdString[ifrq])
            norms1.add_normalMode(norm1)
        vib1.set_normalModes(norms1)
    # dipole moment information
    molDipoleX = psi4.get_variable('CURRENT DIPOLE X')
    molDipoleY = psi4.get_variable('CURRENT DIPOLE Y')
    molDipoleZ = psi4.get_variable('CURRENT DIPOLE Z')
    molDipoleTot = math.sqrt(
        molDipoleX * molDipoleX +
        molDipoleY * molDipoleY +
        molDipoleZ * molDipoleZ)
    prop1 = api.propertiesType()
    sprop1 = api.propertyType(
        name='dipoleMomentX',
        unit='gc:debye')
    sprop1.set_valueOf_(molDipoleX)
    sprop2 = api.propertyType(
        name='dipoleMomentY',
        unit='gc:debye')
    sprop2.set_valueOf_(molDipoleY)
    sprop3 = api.propertyType(
        name='dipoleMomentZ',
        unit='gc:debye')
    sprop3.set_valueOf_(molDipoleZ)
    sprop4 = api.propertyType(
        name='dipoleMomentAverage',
        unit='gc:debye')
    sprop4.set_valueOf_(molDipoleTot)
    prop1.add_systemProperty(sprop1)
    prop1.add_systemProperty(sprop2)
    prop1.add_systemProperty(sprop3)
    prop1.add_systemProperty(sprop4)

    # get the basename for the CSX file
    psio = psi4.IO.shared_object()
    namespace = psio.get_default_namespace()
    #csxfilename = '.'.join([namespace, str(os.getpid()), 'csx'])
    csxfilename = os.path.splitext(psi4.outfile_name())[0] + '.csx'
    csxfile = open(csxfilename, 'w')
    csxVer = psi4.get_global_option('CSX_VERSION')

    # Both CSX versions 0 and 1 depended on the procedures table, which in
    #   turn required the writeCSX function to be in the driver.py file
    #   itself. Starting with 1.5 (1, to run), this dependence is broken and
    #   CSX has been shifted into a plugin.

    # Start to generate CSX elements
    # CSX version 1.5
    if csxVer == 2.0:
        #       import csx1_api as api
        cs1 = api.csType(version='2.0') #5')

        # molPublication section: 1.5
        mp1 = api.mpubType(
            title=psi4.get_global_option('PUBLICATIONTITLE'),
            abstract=psi4.get_global_option('PUBLICATIONABSTRACT'),
            publisher=psi4.get_global_option('PUBLICATIONPUBLISHER'),
            status=['PRELIMINARY', 'DRAFT', 'FINAL'].index(psi4.get_global_option('PUBLICATIONSTATUS')),
            category=psi4.get_global_option('PUBLICATIONCATEGORY'),
            visibility=['PRIVATE', 'PROTECTED', 'PUBLIC'].index(psi4.get_global_option('PUBLICATIONVISIBILITY')),
            tags=psi4.get_global_option('PUBLICATIONTAGS'),
            key=psi4.get_global_option('PUBLICATIONKEY'))
        email = psi4.get_global_option('EMAIL').replace('__', '@')
        mp1.add_author(api.authorType(
            creator=psi4.get_global_option('CORRESPONDINGAUTHOR'),
            type_='gc:CorrespondingAuthor',
            organization=psi4.get_global_option('ORGANIZATION'),
            email=None if email == '' else email))
        #mp1 = api.mpType(
        #   title='', abstract='', publisher='', status=0, category=2, visibility=0, tags='', key='')
        mp1.set_sourcePackage(api.sourcePackageType(name='Psi4', version=psi4.version()))
        #mp1.add_author(api.authorType(creator='', type_='cs:corresponding', organization='', email=''))
        cs1.set_molecularPublication(mp1)

        # molSystem section: 1.5
        ms1 = api.msysType(
            systemCharge=molCharge,
            systemMultiplicity=molMulti, id='s1')
        temp1 = api.dataWithUnitsType(unit='gc:kelvin')
        temp1.set_valueOf_(0.0)  # LAB dispute
        ms1.set_systemTemperature(temp1)
        mol1 = api.moleculeType(id='m1', atomCount=molecule.natom())
        #OBmol1 = api.moleculeType(id='m1', atomCount=atomNum)
        #OBobmol1 = openbabel.OBMol()
        #OBfor iatm in range(atomNum):
        #OB    atomField = atomLine[iatm + 1].split()
        #OB    atmSymbol = atomField[0]
        #OB    xCoord = float(atomField[1])
        #OB    yCoord = float(atomField[2])
        #OB    zCoord = float(atomField[3])
        #OB    obatm = obmol1.NewAtom()
        #OB    obatm.SetAtomicNum(qcdb.periodictable.el2z[atmSymbol.upper()])
        #OB    obatm.SetVector(xCoord, yCoord, zCoord)
        #OBobmol1.ConnectTheDots()
        #OBobmol1.PerceiveBondOrders()
        #OBobmol1.SetTotalSpinMultiplicity(molMulti)
        #OBobmol1.SetTotalCharge(molCharge)
        #OBconv1 = openbabel.OBConversion()
        #OBconv1.SetInAndOutFormats('mol', 'inchi')
        #OBconv1.SetOptions('K', conv1.OUTOPTIONS)
        #OBinchikey = conv1.WriteString(obmol1)
        #OBmol1.set_inchiKey(inchikey.rstrip())
        #OBiatm = 0
        for at in range(molecule.natom()):
            #xCoord1 = api.dataWithUnitsType(unit='cs:angstrom')
            #yCoord1 = api.dataWithUnitsType(unit='cs:angstrom')
            #zCoord1 = api.dataWithUnitsType(unit='cs:angstrom')
            #xCoord1.set_valueOf_(molecule.x(at) * p4const.psi_bohr2angstroms)
            #yCoord1.set_valueOf_(molecule.y(at) * p4const.psi_bohr2angstroms)
            #zCoord1.set_valueOf_(molecule.z(at) * p4const.psi_bohr2angstroms)
            xCoord1 = api.dataWithUnitsType(unit='gc:bohr')
            yCoord1 = api.dataWithUnitsType(unit='gc:bohr')
            zCoord1 = api.dataWithUnitsType(unit='gc:bohr')
            xCoord1.set_valueOf_(molecule.x(at))
            yCoord1.set_valueOf_(molecule.y(at))
            zCoord1.set_valueOf_(molecule.z(at))
            # LAB 8jun2015: not getting masses from OB anymore so now dependent on qc programs
            #   current proposition is changing API so masses only go into CSX if relevant (e.g., vib)
            #   same situation as temperature
            atm = api.atomType(
                id='a' + str(at + 1),
                elementSymbol=molecule.symbol(at),
                atomMass=molecule.mass(at),  # psi4 uses mass of most common isotope; OB uses natural distribution mass
                xCoord3D=xCoord1,
                yCoord3D=yCoord1,
                zCoord3D=zCoord1,
                basisSet='bse:' + molBasis,
                calculatedAtomCharge=0,
                formalAtomCharge=0)
            #OBiatm += 1
            #OBcoord1 = api.coordinationType()
            #OBibond = 0
            #OBfor nb_atom in openbabel.OBAtomAtomIter(obatom):
            #OB    bond = obatom.GetBond(nb_atom)
            #OB    bond1 = api.bondType(
            #OB        id1='a' + str(obatom.GetId() + 1),
            #OB        id2='a' + str(nb_atom.GetId() + 1))
            #OB    if bond.GetBondOrder() == 1:
            #OB        bond1.set_valueOf_('single')
            #OB    elif bond.GetBondOrder() == 2:
            #OB        bond1.set_valueOf_('double')
            #OB    elif bond.GetBondOrder() == 3:
            #OB        bond1.set_valueOf_('triple')
            #OB    elif bond.GetBondOrder() == 5:
            #OB        bond1.set_valueOf_('aromatic')
            #OB    else:
            #OB        print('wrong bond order')
            #OB    coord1.add_bond(bond1)
            #OB    ibond += 1
            #OBcoord1.set_bondCount(ibond)
            #OBatm.set_coordination(coord1)
            mol1.add_atom(atm)
        ms1.add_molecule(mol1)
        cs1.set_molecularSystem(ms1)

        # molCalculation section: 1.5
        mc1 = api.mcalType(id='c1')
        qm1 = api.qmCalcType()
        srs1 = api.srsMethodType()
        psivars = psi4.get_variables()

        def form_ene(mandatoryPsivars, optionalPsivars={}, excessPsivars={}):
            """

            """
            ene = api.energiesType(unit='gc:hartree')
            for pv, csx in mandatoryPsivars.iteritems():
                term = api.energyType(type_=csx)
                term.set_valueOf_(psivars.pop(pv))
                ene.add_energy(term)
            for pv, csx in optionalPsivars.iteritems():
                if pv in psivars:
                    term = api.energyType(type_=csx)
                    term.set_valueOf_(psivars.pop(pv))
                    ene.add_energy(term)
            for pv in excessPsivars:
                if pv in psivars:
                    psivars.pop(pv)
            return ene

        # Reference stage- every calc has one
        if 'CCSD TOTAL ENERGY' in psivars or 'CCSD(T) TOTAL ENERGY' in psivars \
                or 'CISD TOTAL ENERGY' in psivars or 'FCI TOTAL ENERGY' in psivars \
                or 'QCISD TOTAL ENERGY' in psivars or 'QCISD(T) TOTAL ENERGY' in psivars:
            mdm1 = api.srsmdMethodType()
            # CCSD(T): 1.5
            if 'CCSD(T) TOTAL ENERGY' in psivars:
                mandatoryPsivars = {
                    'CCSD(T) CORRELATION ENERGY': 'gc:correlation',
                    'CCSD(T) TOTAL ENERGY': 'gc:electronic'}
                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed CCSD(T) computation""")

                block = api.resultType(  # TODO should be pointing to HF for correlation, maybe to MP2 for guess
                    methodology='gc:normal',  # TODO handle dfcc
                    spinType='gc:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                    basisSet='bse:' + molBasis)
                block.set_energies(form_ene(mandatoryPsivars))
                if hasFreq:
                    block.set_vibrationalAnalysis(vib1)
                mdm1.set_ccsd_t(block)
            # CCSD: 1.5
            elif 'CCSD TOTAL ENERGY' in psivars:
                mandatoryPsivars = {
                    'CCSD CORRELATION ENERGY': 'gc:correlation',
                    'CCSD TOTAL ENERGY': 'gc:electronic'}
                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed CCSD computation""")

                block = api.resultType(  # TODO should be pointing to HF for correlation, maybe to MP2 for guess
                    methodology='gc:normal',  # TODO handle dfcc
                    spinType='gc:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                    basisSet='bse:' + molBasis)
                block.set_energies(form_ene(mandatoryPsivars))
                if hasFreq:
                    block.set_vibrationalAnalysis(vib1)
                mdm1.set_ccsd(block)
            # CISD: 1.5
            elif 'CISD TOTAL ENERGY' in psivars:
                mandatoryPsivars = {
                    'CISD CORRELATION ENERGY': 'gc:correlation',
                    'CISD TOTAL ENERGY': 'gc:electronic'}
                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed CISD computation""")

                block = api.resultType(  # TODO should be pointing to HF for correlation, maybe to MP2 for guess
                    methodology='gc:normal',  # TODO handle dfcc
                    spinType='gc:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                    basisSet='bse:' + molBasis)
                block.set_energies(form_ene(mandatoryPsivars))
                if hasFreq:
                    block.set_vibrationalAnalysis(vib1)
                mdm1.set_cisd(block)
            # FCI: 1.5
            elif 'FCI TOTAL ENERGY' in psivars:
                mandatoryPsivars = {
                    'FCI CORRELATION ENERGY': 'gc:correlation',
                    'FCI TOTAL ENERGY': 'gc:electronic'}
                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed FCI computation""")

                block = api.resultType(  # TODO should be pointing to HF for correlation, maybe to MP2 for guess
                    methodology='gc:normal',  # TODO handle dfcc
                    spinType='gc:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                    basisSet='bse:' + molBasis)
                block.set_energies(form_ene(mandatoryPsivars))
                mdm1.set_fci(block)
            # QCISD(T): 1.5
            elif 'QCISD(T) TOTAL ENERGY' in psivars:
                mandatoryPsivars = {
                    'QCISD(T) CORRELATION ENERGY': 'gc:correlation',
                    'QCISD(T) TOTAL ENERGY': 'gc:electronic'}
                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed QCISD(T) computation""")

                block = api.resultType(  # TODO should be pointing to HF for correlation, maybe to MP2 for guess
                    methodology='gc:normal',  # TODO handle dfcc
                    spinType='gc:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                    basisSet='bse:' + molBasis)
                block.set_energies(form_ene(mandatoryPsivars))
                if hasFreq:
                    block.set_vibrationalAnalysis(vib1)
                mdm1.set_qcisd_t(block)
            # QCISD: 1.5
            elif 'QCISD TOTAL ENERGY' in psivars:
                mandatoryPsivars = {
                    'QCISD CORRELATION ENERGY': 'gc:correlation',
                    'QCISD TOTAL ENERGY': 'gc:electronic'}
                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed QCISD computation""")

                block = api.resultType(  # TODO should be pointing to HF for correlation, maybe to MP2 for guess
                    methodology='gc:normal',  # TODO handle dfcc
                    spinType='gc:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                    basisSet='bse:' + molBasis)
                block.set_energies(form_ene(mandatoryPsivars))
                if hasFreq:
                    block.set_vibrationalAnalysis(vib1)
                mdm1.set_qcisd(block)
            srs1.set_multipleDeterminant(mdm1)

        elif 'DFT TOTAL ENERGY' in psivars or 'HF TOTAL ENERGY' in psivars \
                or 'MP2 TOTAL ENERGY' in psivars or 'MP3 TOTAL ENERGY' in psivars \
                or 'MP4 TOTAL ENERGY' in psivars:
            sdm1 = api.srssdMethodType()
            # DFT 1.5
            if 'DFT TOTAL ENERGY' in psivars:  # TODO robust enough to avoid MP2C, etc.?
                mandatoryPsivars = {
                    'NUCLEAR REPULSION ENERGY': 'gc:nuclearRepulsion',
                    'DFT FUNCTIONAL TOTAL ENERGY': 'gc:dftFunctional',
                    'DFT TOTAL ENERGY': 'gc:electronic'}
                optionalPsivars = {
                    'DOUBLE-HYBRID CORRECTION ENERGY': 'gc:doubleHybrid correction',
                    'DISPERSION CORRECTION ENERGY': 'gc:dispersion correction'}
                excessPsivars = [
                    'MP2 TOTAL ENERGY',
                    'MP2 CORRELATION ENERGY',
                    'MP2 SAME-SPIN CORRELATION ENERGY']

                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed DFT computation""")

                block = api.resultType(
                    methodology='gc:normal',  # TODO handle dfhf, dfmp
                    spinType='gc:' + molSpin,
                    basisSet='bse:' + molBasis,
                    dftFunctional=name)  # TODO this'll need to be exported
                block.set_energies(form_ene(mandatoryPsivars, optionalPsivars, excessPsivars))
                if wfn:
                    block.set_waveFunction(wfn1)
                if hasFreq:
                    block.set_vibrationalAnalysis(vib1)
                block.set_properties(prop1)
                sdm1.set_dft(block)


            # post-reference block
            # MP4: 1.5
            elif 'MP4 TOTAL ENERGY' in psivars:
                mandatoryPsivars = {
                    'MP4 CORRELATION ENERGY': 'gc:correlation',
                    'MP4 TOTAL ENERGY': 'gc:electronic'}
                optionalPsivars = {
                    'MP4 SAME-SPIN CORRELATION ENERGY': 'gc:sameSpin correlation'}
                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed MP4 computation""")

                block = api.resultType(  # TODO should be pointing to HF for correlation
                    methodology='gc:normal',  # TODO handle dfmp
                    spinType='gc:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                    basisSet='bse:' + molBasis)
                block.set_energies(form_ene(mandatoryPsivars, optionalPsivars))
                if wfn:
                    block.set_waveFunction(wfn1)
                if hasFreq:
                    block.set_vibrationalAnalysis(vib1)
                block.set_properties(prop1)
                sdm1.set_mp4(block)

            # MP3: 1.5
            elif 'MP3 TOTAL ENERGY' in psivars:
                mandatoryPsivars = {
                    'MP3 CORRELATION ENERGY': 'gc:correlation',
                    'MP3 TOTAL ENERGY': 'gc:electronic'}
                optionalPsivars = {
                    'MP3 SAME-SPIN CORRELATION ENERGY': 'gc:sameSpin correlation'}
                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed MP3 computation""")

                block = api.resultType(  # TODO should be pointing to HF for correlation
                    methodology='gc:normal',  # TODO handle dfmp
                    spinType='gc:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                    basisSet='bse:' + molBasis)
                block.set_energies(form_ene(mandatoryPsivars, optionalPsivars))
                if wfn:
                    block.set_waveFunction(wfn1)
                if hasFreq:
                    block.set_vibrationalAnalysis(vib1)
                block.set_properties(prop1)
                sdm1.set_mp3(block)

            # MP2: 1.5
            elif 'MP2 TOTAL ENERGY' in psivars:
                mandatoryPsivars = {
                    'MP2 CORRELATION ENERGY': 'gc:correlation',
                    'MP2 TOTAL ENERGY': 'gc:electronic'}
                optionalPsivars = {
                    'MP2 SAME-SPIN CORRELATION ENERGY': 'gc:sameSpin correlation'}
                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed MP2 computation""")

                block = api.resultType(  # TODO should be pointing to HF for correlation
                    methodology='gc:normal',  # TODO handle dfmp
                    spinType='gc:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                    basisSet='bse:' + molBasis)
                block.set_energies(form_ene(mandatoryPsivars, optionalPsivars))
                if wfn:
                    block.set_waveFunction(wfn1)
                if hasFreq:
                    block.set_vibrationalAnalysis(vib1)
                block.set_properties(prop1)
                sdm1.set_mp2(block)

            # SCF: 1.5
            elif 'HF TOTAL ENERGY' in psivars:
                mandatoryPsivars = {
                    'NUCLEAR REPULSION ENERGY': 'gc:nuclearRepulsion',
                    'HF TOTAL ENERGY': 'gc:electronic'}

                if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                    raise CSXError("""Malformed HF computation""")

                block = api.resultType(
                    methodology='gc:normal',  # TODO handle dfhf, dfmp
                    spinType='gc:' + molSpin,
                    basisSet='bse:' + molBasis)
                block.set_energies(form_ene(mandatoryPsivars))
                if wfn:
                    block.set_waveFunction(wfn1)
                if hasFreq:
                    block.set_vibrationalAnalysis(vib1)
                block.set_properties(prop1)
                sdm1.set_abinitioScf(block)
            else:
                psi4.print_out("""\nCSX version {0} does not support """
                               """method {1} for {2}\n""".format(
                               csxVer, lowername, 'energies'))

            srs1.set_singleDeterminant(sdm1)

        #print('CSX not harvesting: ', ', '.join(psivars))

        qm1.set_singleReferenceState(srs1)
        mc1.set_quantumMechanics(qm1)
        cs1.set_molecularCalculation(mc1)

    else:
        print('The future CSX file is here')

    csxfile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    cs1.export(csxfile, 0)
    csxfile.close()
    # End to write the CSX file

import procedures
procedures.proc_table.hooks['energy']['post'].append(writeCSX)
procedures.proc_table.hooks['optimize']['post'].append(writeCSX)
procedures.proc_table.hooks['frequency']['post'].append(writeCSX)

