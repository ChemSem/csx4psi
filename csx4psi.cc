/*
 *@BEGIN LICENSE
 *
 * csx4psi by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace csx4psi {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "CSX4PSI"|| options.read_globals()) {
        /*- Do write a CSX output file? -*/
        options.add_bool("WRITE_CSX", true);
        /*- Version of CSX. -*/
        options.add_int("CSX_VERSION", 1);
        /*- Brief title for CSX publication. -*/
        options.add_str("PUBLICATIONTITLE", "");
//        /*- Abstract for CSX publication. -*/
        options.add_str("PUBLICATIONABSTRACT", "");
//        /*- Additional tags to aid publication searches. -*/
        options.add_str("PUBLICATIONTAGS", "");
//        /*- Submitting author for CSX publication. -*/
        options.add_str("publicationPublisher", "");
//        /*- Corresponding author for CSX publication. -*/
        options.add_str("CORRESPONDINGAUTHOR", "");
//        /*- Organization to which corresponding or submitting author belongs. -*/
        options.add_str("ORGANIZATION", "");
//        /*- E-mail address of record for CSX publication -*/
        options.add_str("EMAIL", "");
//        /*- Area of chemistry. -*/
        options.add_int("PUBLICATIONCATEGORY", 2);
//        /*- Visibility of CSX publication. ``PRIVATE`` accessible only by
//        submitting author with password. ``PROTECTED`` accessible with key
//        provided by submitting author. ``PUBLIC`` freely accessible. -*/
        options.add_str("PUBLICATIONVISIBILITY", "PRIVATE", "PRIVATE PROTECTED PUBLIC");
//        /*- Status of CSX publication. ``PRELIMINARY`` is unreviewed.
//        ``DRAFT`` is not fully reviewed. ``FINAL`` is final. -*/
        options.add_str("PUBLICATIONSTATUS", "PRELIMINARY", "PRELIMINARY DRAFT FINAL");
//        /*-  -*/
        options.add_int("PUBLICATIONKEY", 0);
    }

    return true;
}

extern "C" 
SharedWavefunction csx4psi(SharedWavefunction ref_wfn, Options& options)
{
    int print = options.get_int("PRINT");

    /* Your code goes here */

    return ref_wfn;
}

}} // End namespaces

