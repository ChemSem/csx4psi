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
from driver import *
from wrappers import *
from molutil import *
import p4util
from p4xcpt import *


def writeCSX(name, **kwargs):
    """function to write the CSX file

    """

    if not psi4.get_global_option('WRITE_CSX'):
        return

    # import csx_api for csx writing
    import os
    import math
    import inspect
    import openbabel
    import qcdb
    import qcdb.periodictable
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
        'hessian': 2,
        }
    dertype = derdict[calledby]
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
        wfn = psi4.wavefunction()
    except AttributeError:
        pass
    if wfn:
        molOrbE = psi4.wavefunction().epsilon_a()
        molOrbEb = psi4.wavefunction().epsilon_b()
        orbNmopi = psi4.wavefunction().nmopi()
        orbNsopi = psi4.wavefunction().nsopi()
        orbNum = psi4.wavefunction().nmo() if molOrbE else 0
        orbSNum = psi4.wavefunction().nso()
        molOrb = psi4.wavefunction().Ca()
        orbNirrep = psi4.wavefunction().nirrep()
        orbAotoso = psi4.wavefunction().aotoso()
        orbDoccpi = psi4.wavefunction().doccpi()
        orbSoccpi = psi4.wavefunction().soccpi()
        basisNbf = psi4.wavefunction().basisset().nbf()
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
            molOrbCb = psi4.wavefunction().Cb()
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
    # frequency information
    if dertype == 2:
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
    molDipoleX = psi4.get_variable('CURRENT DIPOLE X')
    molDipoleY = psi4.get_variable('CURRENT DIPOLE Y')
    molDipoleZ = psi4.get_variable('CURRENT DIPOLE Z')
    molDipoleTot = math.sqrt(
        molDipoleX * molDipoleX +
        molDipoleY * molDipoleY +
        molDipoleZ * molDipoleZ)

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
    if csxVer == 1.0:
        import csx1_api as api
        cs1 = api.csType(version='1.0') #5')

        # molPublication section: 1.5
        mp1 = api.mpType(
            title=psi4.get_global_option('PUBLICATIONTITLE'),
            abstract=psi4.get_global_option('PUBLICATIONABSTRACT'),
            publisher=psi4.get_global_option('PUBLICATIONPUBLISHER'),
            status=['PRELIMINARY', 'DRAFT', 'FINAL'].index(psi4.get_global_option('PUBLICATIONSTATUS')),
            category=psi4.get_global_option('PUBLICATIONCATEGORY'),
            visibility=['PRIVATE', 'PROTECTED', 'PUBLIC'].index(psi4.get_global_option('PUBLICATIONVISIBILITY')),
            tags=psi4.get_global_option('PUBLICATIONTAGS'),
            key=psi4.get_global_option('PUBLICATIONKEY'))
        mp1.set_sourcePackage(api.sourcePackageType(name='Psi4', version=psi4.version()))
        email = psi4.get_global_option('EMAIL').replace('__', '@')
        mp1.add_author(api.authorType(
            creator=psi4.get_global_option('CORRESPONDINGAUTHOR'),
            type_='cs:corresponding',
            organization=psi4.get_global_option('ORGANIZATION'),
            email=None if email == '' else email))
        cs1.set_molecularPublication(mp1)

        # molSystem section: 1.5
        ms1 = api.msType(
            systemCharge=molCharge,
            systemMultiplicity=molMulti)
        temp1 = api.dataWithUnitsType(unit='cs:kelvin')
        temp1.set_valueOf_(298.0)  # LAB dispute
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
            xCoord1 = api.dataWithUnitsType(unit='cs:bohr')
            yCoord1 = api.dataWithUnitsType(unit='cs:bohr')
            zCoord1 = api.dataWithUnitsType(unit='cs:bohr')
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
                basisSet='cs:' + molBasis,
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
        avalMethods = False
        mc1 = api.mcType()
        qm1 = api.qmCalcType()
        srs1 = api.srsMethodType()
        sdm1 = api.srssdMethodType()
        psivars = psi4.get_variables()

        def form_ene(mandatoryPsivars, optionalPsivars={}, excessPsivars={}):
            """

            """
            ene = api.energiesType(unit='cs:hartree')
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
        # DFT 1.5
        if 'DFT TOTAL ENERGY' in psivars:  # TODO robust enough to avoid MP2C, etc.?
            mandatoryPsivars = {
                'NUCLEAR REPULSION ENERGY': 'cs:nuclearRepulsion',
                'DFT FUNCTIONAL TOTAL ENERGY': 'cs:dftFunctional',
                'DFT TOTAL ENERGY': 'cs:electronic'}
            optionalPsivars = {
                'DOUBLE-HYBRID CORRECTION ENERGY': 'cs:doubleHybrid correction',
                'DISPERSION CORRECTION ENERGY': 'cs:dispersion correction'}
            excessPsivars = [
                'MP2 TOTAL ENERGY',
                'MP2 CORRELATION ENERGY',
                'MP2 SAME-SPIN CORRELATION ENERGY']

            if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                raise CSXError("""Malformed DFT computation""")

            block = api.resultType(
                methodology='cs:normal',  # TODO handle dfhf, dfmp
                spinType='cs:' + molSpin,
                basisSet='bse:' + molBasis,
                dftFunctional=name)  # TODO this'll need to be exported
            block.set_energies(form_ene(mandatoryPsivars, optionalPsivars, excessPsivars))
            sdm1.set_dft(block)

        # SCF: 1.5
        elif 'HF TOTAL ENERGY' in psivars:
            mandatoryPsivars = {
                'NUCLEAR REPULSION ENERGY': 'cs:nuclearRepulsion',
                'HF TOTAL ENERGY': 'cs:electronic'}

            if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                raise CSXError("""Malformed HF computation""")

            block = api.resultType(
                methodology='cs:normal',  # TODO handle dfhf, dfmp
                spinType='cs:' + molSpin,
                basisSet='bse:' + molBasis)
            block.set_energies(form_ene(mandatoryPsivars))
            sdm1.set_abinitioScf(block)
        else:
            psi4.print_out("""\nCSX version {0} does not support """
                           """method {1} for {2}\n""".format(
                           csxVer, lowername, 'energies'))

        # post-reference block
        # MP2: 1.5
        if 'MP2 TOTAL ENERGY' in psivars:
            mandatoryPsivars = {
                'MP2 CORRELATION ENERGY': 'cs:correlation'}
            optionalPsivars = {
                'MP2 SAME-SPIN CORRELATION ENERGY': 'cs:sameSpin correlation'}
            if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                raise CSXError("""Malformed MP2 computation""")

            block = api.resultType(  # TODO should be pointing to HF for correlation
                methodology='cs:normal',  # TODO handle dfmp
                spinType='cs:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                basisSet='bse:' + molBasis)
            block.set_energies(form_ene(mandatoryPsivars, optionalPsivars))
            sdm1.set_mp2(block)

        # CCSD: 1.5
        if 'CCSD TOTAL ENERGY' in psivars:
            mandatoryPsivars = {
                'CCSD CORRELATION ENERGY': 'cs:correlation'}
            optionalPsivars = {
                'CCSD SAME-SPIN CORRELATION ENERGY': 'cs:sameSpin correlation'}
            excessPsivars = [
                'CCSD TOTAL ENERGY',
                'CCSD OPPOSITE-SPIN CORRELATION ENERGY']
            if not all([pv in psivars for pv in mandatoryPsivars.keys()]):
                raise CSXError("""Malformed CCSD computation""")

            block = api.resultType(  # TODO should be pointing to HF for correlation, maybe to MP2 for guess
                methodology='cs:normal',  # TODO handle dfcc
                spinType='cs:' + molSpin,  # TODO could have a closed-shell corl mtd atop open-shell scf?
                basisSet='bse:' + molBasis)
            block.set_energies(form_ene(mandatoryPsivars, optionalPsivars))
            #sdm1.set_ccsd(block)  # needs to be srsmd

        #print('CSX not harvesting: ', ', '.join(psivars))
        # LAB TODO not addressed below here
        # wavefunction: 1.5
        if avalMethods:
            if wfnRestricted:
                wfn1 = api.waveFunctionType(
                    orbitalCount=orbNum,
                    orbitalOccupancies=orbOccString)
                orbe1 = api.stringArrayType(unit='cs:hartree')
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
                orbe1 = api.stringArrayType(unit='cs:hartree')
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
                orbeb1 = api.stringArrayType(unit='cs:hartree')
                orbeb1.set_valueOf_(orbEbString)
                wfn1.set_betaOrbitalEnergies(orbeb1)
                wfn1.set_betaOrbitalOccupancies(orbOccCbString)
                borbs1 = api.orbitalsType()
                for iorb in range(orbNum):
                    orb1 = api.stringArrayType(id=iorb+1)
                    orb1.set_valueOf_(orbCbString[iorb])
                    borbs1.add_orbital(orb1)
                wfn1.set_betaOrbitals(borbs1)

            scf1.set_waveFunction(wfn1)
            if dertype == 2:
                vib1 = api.vibAnalysisType(vibrationCount=molFreqNum)
                freq1 = api.stringArrayType(unit="cs:cm-1")
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
                scf1.set_vibrationalAnalysis(vib1)
            # Properties: 1.5
            prop1 = api.propertiesType()
            sprop1 = api.propertyType(
                name='dipoleMomentX',
                unit='cs:debye')
            sprop1.set_valueOf_(molDipoleX)
            sprop2 = api.propertyType(
                name='dipoleMomentY',
                unit='cs:debye')
            sprop2.set_valueOf_(molDipoleY)
            sprop3 = api.propertyType(
                name='dipoleMomentZ',
                unit='cs:debye')
            sprop3.set_valueOf_(molDipoleZ)
            sprop4 = api.propertyType(
                name='dipoleMomentAverage',
                unit='cs:debye')
            sprop4.set_valueOf_(molDipoleTot)
            prop1.add_systemProperty(sprop1)
            prop1.add_systemProperty(sprop2)
            prop1.add_systemProperty(sprop3)
            prop1.add_systemProperty(sprop4)
            scf1.set_properties(prop1)

        srs1.set_singleDeterminant(sdm1)
        qm1.set_singleReferenceState(srs1)
        mc1.set_quantumMechanics(qm1)
        cs1.set_molecularCalculation(mc1)

    else:
        print('The future CSX file is here')

    csxfile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    cs1.export(csxfile, 0)
    csxfile.close()
    # End to write the CSX file

hooks['energy']['post'].append(writeCSX)

