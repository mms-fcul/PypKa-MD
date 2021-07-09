/*
This file is part of ST-CpHMD, version v4.1_GMX2018.

Copyright (c) 2005-2020, Instituto de Tecnologia Quimica e Biologica,
Universidade Nova de Lisboa, Portugal.

ST-CpHMD is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 2 of the License, or (at your
option) any later version.

ST-CpHMD is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with ST-CpHMD.  If not, see <http://www.gnu.org/licenses/>.

For further details and info check the manual.pdf file.

You can get ST-CpHMD at www.itqb.unl.pt/simulation
*/


/*

                     *** FixBox, version 1.2 ***

  This program reads a .gro file and, for each frame, assembles and centers
  the system in the periodic box, using information provided in an input
  file with molecular definitions (moldef), writing the result to stdout.

  The moldef file is read only after the first frame and contains the
  definition of all the necessary groups of molecules (which together must
  define *all* molecules), the indication of the groups to be used for
  assembling and centering, and the indication of the final PBCs (see
  format below). The part containing the group definitions is a
  quick-and-dirty substitute for a molecular topology and its parsing is
  somewhat cumbersome (because of the strict variable declaration and
  allocation in C), but it gets the job done.

  For each frame, the following sequential algorithm is applied:

  0. Determine distances between molecules, using their closest atom.

  1. Assemble the groups of interest by molecular proximity, gradually
     adding the closest molecules, one at a time. More than one stage can
     be used, for greater flexibility.

  2. Center the groups of interest along each direction, using the total
     extent.

  3. Apply boundary conditions to specific groups along each direction.

  Boxes are always assumed to be triclinic. However, except for step 0,
  which is done in physical space, all other steps are performed in the
  corresponding dimensionless unit cubic box in scaled space [Tuckerman
  (2010) "Statistical Mechanics", secs. 5.6 and 5.7, and Appendix B].

  ----------------------------------------------------------------------
  MOLDEF FORMAT:

  The processed lines of the moldef file must start with one of the
  following characters (which usually occur in the file in this order):

  - G : A group name definition, consisting of a string (without
    spaces). Line examples: "G Protein", "G Lipids". The group 'None' is
    predefined and reserved, corresponding to a group with no molecules
    that can be used in other definitions (see below). A group definition
    is followed by the definition of its molecules, one per line, starting
    with one of the lower-case characters listed below.

  - a : A molecule defined by the _ordinal_ range of its atoms. Line
    example: "a 3521 5283", which means that a molecule starts at the
    3521st atom and ends at the 5283rd one, regardless of their numbers in
    the .gro file. The molecule is assigned to the currently defined group
    (as defined by the last G-line).

  - n : One or more molecules defined by the residue name. Line examples:
    "n ATP", "n DMPC". The molecules are assigned to the currently defined
    group.

  - g : One or more molecules defined through a previously defined group
    (not the current one). Line examples: "g Protein", "g Lipids". This is
    useful to define larger groups (e.g., a group System containing the
    groups Protein and Lipids). The molecules are assigned to the currently
    defined group.

  - A : Specify a group to be used (sequentially) in an assembling stage,
    giving one line per stage (1st A-line for stage 1, 2nd A-line for stage
    2, etc). Line examples: "A Protein", "A Lipids".

  - C : Specify the three groups to be centered along each box vector
    (a,b,c), and the corresponding types of message, 'E' or 'W' (error or
    warning), to generate if the group exceeds the box along that
    direction. Line example: "C Protein Protein System E E W", which
    centers Protein along directions a and b, and System along direction c,
    giving an error if Protein exceeds the extent of a and/or b but
    allowing System to exceed the c extent (writing just a
    warning). Together with the P-line, this E/W selection provides
    detailed control over the centering/PBC behavior, avoiding to
    accidentally apply PBCs to groups intended to remain whole.

  - P : Specify the three groups for which PBCs should be applied along
    each box vector (a,b,c). Line example: "P System System None", which
    means applying PBCs to group System along directions a and b, but do
    not apply PBCs to any group along direction c.

  - Lines not starting with one of the above characters are ignored.

  For clarity in the moldef file, it is suggested that: upper-case entries
  are given in the above order (G,A,C,P); upper-case entries are separated
  using blank lines, except A-lines (which should be consecutive);
  lower-case entries follow immediately after their corresponding G
  definition, using the above order (a,n,g); lines with comments are
  started using a character commonly used for that purpose (e.g., '#'). See
  the moldef files provided as examples.

  ----------------------------------------------------------------------
  IMPORTANT CONDITIONS:

  - The input coordinates MUST correspond to whole molecules (e.g.,
    generated with '-pbc whole' or '-pbc mol' by GROMACS). This is NOT
    (and, without a topology, cannot be) checked!

  - A molecule MUST be a set of consecutive atoms. This is NOT checked!

  - As noted above, a-lines in the moldef file should contain _ordinal_
    atom numbers. Actual atom numbers in the .gro file are treated as mere
    string labels and reused only to output the fixed .gro file.

  - The moldef definitions should assign all atoms to a molecule and all
    molecules to at least one group. Otherwise, the program reports an
    error and exits.

  ----------------------------------------------------------------------
  COMPILATION (examples):

    gcc fixbox.c -o fixbox -O3 -lm -Wall -W -pedantic -Wno-unused-result
    gcc fixbox.c -o fixbox -O3 -lm -Wall -W -pedantic -ansi

  ----------------------------------------------------------------------
  AUTHORS:

  - Antonio M. Baptista, ITQB NOVA, Portugal

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <math.h>

#define MAXSTR 200
#define BIG_DISTANCE 1e20
#define TINY_VOLUME 1e-10

typedef struct {
  int anumb, rnumb, mol ;
  float r[3], s[3] ;
  char aname[6+1], rname[6+1], vels[24+1] ;
} atom_t ;

typedef struct {
  int nmols, *mol ;
  char name[MAXSTR] ;
} group_t ;

typedef struct {
  int firstatom, lastatom ;
} mol_t ;

typedef struct {
  float dRMI2, dS[3] ;
} dist_t ;

typedef struct {
  int mol, dN[3], assembled ;
} obj_t ;


/* Function prototypes */
void parse_arguments(int argc, char **argv) ;
void read_frame(void) ;
void read_moldef(void) ;
void write_frame(void) ;
int gnumb(char groupname[]) ;
void r_to_s(void) ;
void s_to_r(void) ;
void compute_objdist(void) ;
void assemble_objs(void) ;
void center(void) ;
void apply_PBCs(void) ;
int nint(float x) ;
void message(char mtype, char *format, ...) ;

/* Global variables */
int natoms, ngroups, nmols, nstages, *nobjs, nobjs_tot, gcent[3],
  gPBC[3], nframe = 0, nClines = 0, nPlines = 0 ;
float H[3][3], iH[3][3], boxC[3] ;  // maybe collect into a box structure
char title[MAXSTR+1], boxline[MAXSTR+1], centEW[3] ;
static const char empty_gname[] = "None" ;
atom_t *atom ;
group_t *group ;
mol_t *mol ;
dist_t **dist ;
obj_t *obj ;
FILE *fp_gro, *fp_moldef ;


int main(int argc, char **argv)
{
  parse_arguments(argc, argv) ;

  while (fgets(title, sizeof title, fp_gro) != NULL)
  {
    nframe++ ;
    read_frame() ;
    if (nframe == 1) read_moldef() ;
    r_to_s() ;
    compute_objdist() ;
    assemble_objs() ;
    center() ;
    apply_PBCs() ;
    s_to_r() ;
    write_frame() ;
  }
  fclose(fp_gro) ;

  // free stuff...

  return 0 ;
}


/*Reads arguments*/
void parse_arguments(int argc, char **argv)
{
  char *gro_file, *moldef_file ;
  
  if (argc != 3) message('U', "Wrong number of arguments.\n") ;

  gro_file = calloc(strlen(argv[1])+1, sizeof(char)) ;
  strcpy(gro_file, argv[1]) ;
  if ((fp_gro = fopen(gro_file, "r")) == NULL)
    message('E', "Opening \"%s\": %s\n", gro_file, strerror(errno));

  moldef_file = calloc(strlen(argv[2])+1, sizeof(char)) ;
  strcpy(moldef_file, argv[2]) ;
  if ((fp_moldef = fopen(moldef_file, "r")) == NULL)
    message('E', "Opening \"%s\": %s\n", moldef_file, strerror(errno));
}


void read_frame(void)
{
  int a, k, j ;
  float detH ;
  char line[MAXSTR+1], resname[6], atmname[6] ;

  //fgets(title, sizeof title, fp_gro) ;

  fgets(line, sizeof line, fp_gro) ;
  sscanf(line, "%d", &natoms) ;
  atom = calloc(natoms+1, sizeof(atom_t)) ;

  // Read atoms:
  for (a = 1 ; a <= natoms ; a++)
  {
    fgets(line, sizeof line, fp_gro) ;
    sscanf(line, "%5d%5[^\n]%5[^\n]%5d%8f%8f%8f%24[^\n]",
	   &(atom[a].rnumb), resname, atmname, &(atom[a].anumb),
	   &(atom[a].r[0]), &(atom[a].r[1]), &(atom[a].r[2]), atom[a].vels) ;
    sscanf(resname, "%s", atom[a].rname) ;
    sscanf(atmname, "%s", atom[a].aname) ;
  }
  fgets(boxline, sizeof boxline, fp_gro) ;
  // Read box matrix H (free format space-separated reals, with 3 or 9 reals):
  // Order = a(x) b(y) c(z) a(y) a(z) b(x) b(z) c(x) c(y)
  // If the last 6 fields do not exist, they should be set to zero.
  // The usual actual format in GROMACS seems to be " %9.5" for each.
  sscanf(boxline, "%f %f %f %f %f %f %f %f %f",
	 &(H[0][0]), &(H[1][1]), &(H[2][2]), &(H[1][0]), &(H[2][0]),
	 &(H[0][1]), &(H[2][1]), &(H[0][2]), &(H[1][2])) ;

  // Compute inverse of H using the determinant+cofactors method:
  detH = H[0][0] * (H[1][1] * H[2][2] - H[1][2] * H[2][1])
         - H[0][1] * (H[1][0] * H[2][2] - H[1][2] * H[2][0])
         + H[0][2] * (H[1][0] * H[2][1] - H[1][1] * H[2][0]) ;
  if (fabs(detH) < TINY_VOLUME)
    message('E', "Tiny box in frame %d. Box info may be missing...\n", nframe);
  for (k = 0 ; k < 3 ; k++)
    for (j = 0 ; j < 3 ; j++)
      iH[k][j] = ( H[(j+1)%3][(k+1)%3] * H[(j+2)%3][(k+2)%3]
	           - H[(j+1)%3][(k+2)%3] * H[(j+2)%3][(k+1)%3] ) / detH ;

  // Compute box center, assuming that the origin is at the "lowest" box corner:
  for (k = 0 ; k < 3 ; k++)
    boxC[k] = (H[k][0] + H[k][1] + H[k][2]) / 2 ;
}


// Strictly speaking, all calloc() calls should be followed by
// assignment of its variables, because zero-bytes do not necessarily
// correspond to zero integers or floats, or even to NULL pointers;
// thus, calloc() is really not much different from malloc() in this
// respect. However, this is perhaps a bit paranoid, since virtually
// all existing compilers actually do that...

// Consider using fseek() instead of rewind().

void read_moldef(void)
{
  // Try to simplify/reduce the number of auxiliary variables...
  int k, g, a, m, gm, prevrnumb, i, subg, st, gst, notot, o ;
  char line[MAXSTR+1], auxname[MAXSTR], auxstr[3][MAXSTR] ;

  // 1st pass: count groups and assembling stages
  ngroups = nstages = 0 ;
  while(fgets(line, sizeof line, fp_moldef) != NULL)
  {
    if (line[0] == 'G') ngroups++ ;
    else if (line[0] == 'A') nstages++ ;
  }
  group = calloc(ngroups+1, sizeof(group_t)) ; // indices = {0 to ngroups}
  nobjs = calloc(nstages, sizeof(int)) ;  // indices = {0 to nstages-1}
  // Define empty group:
  strcpy(group[0].name, empty_gname) ;
  group[0].nmols = 0 ;
  
  // 2nd pass: read group names, count total and per-group molecules
  rewind(fp_moldef) ;
  g = 0 ;
  nmols = 0 ;
  prevrnumb = 0 ;
  st = 0 ;
  nobjs_tot = 0 ;
  while(fgets(line, sizeof line, fp_moldef) != NULL)
  {
    if (line[0] == 'G')
    {
      sscanf(line, "G %s", group[++g].name) ;
      if (strcmp(group[g].name, group[0].name) == 0)
	message('E', "Name '%s' not allowed for user-defined groups.\n",
		group[0].name) ;;
    }
    else if (line[0] == 'a')
    {
      if (g == 0) message('E', "Entry 'a' found before a group was defined.\n") ;
      group[g].nmols++ ;
      nmols++ ;
    }
    else if (line[0] == 'n')
    {
      if (g == 0) message('E', "Entry 'n' found before a group was defined.\n") ;
      sscanf(line, "n %s", auxname) ;
      for (a = 1 ; a <= natoms ; a++)
      {
	if (strcmp(auxname, atom[a].rname) == 0)
	{
	  if (atom[a].rnumb != prevrnumb)
	  {
	    group[g].nmols++ ;
	    nmols++ ;
	  }
	}
	prevrnumb = atom[a].rnumb ;
      }
    }
    else if (line[0] == 'g')
    {
      if (g == 0) message('E', "Entry 'g' found before a group was defined.\n") ;
      sscanf(line, "g %s", auxname) ;
      group[g].nmols += group[gnumb(auxname)].nmols ;
    }
    else if (line[0] == 'A')
    {
      sscanf(line, "A %s", auxname) ;
      nobjs[++st] = nobjs_tot += group[gnumb(auxname)].nmols ;
    }
  }
  mol = calloc(nmols+1, sizeof(mol_t)) ;
  for (g = 1 ; g <= ngroups ; g++)
    group[g].mol = calloc(group[g].nmols+1, sizeof(int)) ;
  obj = calloc(nobjs_tot+1, sizeof(obj_t)) ;
  dist = calloc(nobjs_tot+1, sizeof(dist_t *)) ;
  for (o = 1 ; o <= nobjs_tot ; o++)
    dist[o] = calloc(nobjs_tot+1, sizeof(dist_t)) ;

  // 3rd pass: assign molecules to groups and atoms to molecules.
  rewind(fp_moldef) ;
  g = m = gm = 0 ;
  prevrnumb = 0 ;
  notot = 0 ;
  while(fgets(line, sizeof line, fp_moldef) != NULL)
  {
    if (line[0] == 'G')
    {
      g++ ;
      gm = 0 ;
    }
    else if (line[0] == 'a')
    {
      group[g].mol[++gm] = ++m ;
      sscanf(line, "a %d %d", &(mol[m].firstatom), &(mol[m].lastatom)) ;
      for (a = mol[m].firstatom ; a <= mol[m].lastatom ; a++)
      {
	if (atom[a].mol == 0) atom[a].mol = m ;
	else message('E', "Atom index %d already assigned to molecule.\n", a) ;
      }
    }
    else if (line[0] == 'n')
    {
      sscanf(line, "n %s", auxname) ;
      for (a = 1 ; a <= natoms ; a++)
      {
	if (strcmp(auxname, atom[a].rname) == 0)
	{
	  if (atom[a].rnumb != prevrnumb) group[g].mol[++gm] = ++m ;
	  if (mol[m].firstatom == 0) mol[m].firstatom = a ;
	  mol[m].lastatom = a ;
	  if (atom[a].mol == 0) atom[a].mol = m ;
	  else message('E', "Atom %d already assigned to molecule.\n",
		       atom[a].anumb) ;
	}
	prevrnumb = atom[a].rnumb ;
      }
    }
    else if (line[0] == 'g')
    {
      sscanf(line, "g %s", auxname) ;
      subg = gnumb(auxname) ;
      for (i = 1 ; i <= group[subg].nmols ; i++)
	group[g].mol[++gm] = group[subg].mol[i] ;
    }
    else if (line[0] == 'A')
    {
      sscanf(line, "A %s", auxname) ;
      gst = gnumb(auxname) ;
      for (i = 1 ; i <= group[gst].nmols ; i++)
	obj[++notot].mol = group[gst].mol[i] ;
    }
    else if (line[0] == 'C')
    {
      nClines++ ;
      sscanf(line, "C %s %s %s %c %c %c", auxstr[0], auxstr[1], auxstr[2],
	     &(centEW[0]), &(centEW[1]), &(centEW[2])) ;
      for (k = 0 ; k < 3 ; k++)
      {
	gcent[k] = gnumb(auxstr[k]) ;
	if (centEW[k] != 'E' && centEW[k] != 'W')
	  message('E', "Give either 'E' or 'W' as message codes in C-line.\n") ;
      }
    }
    else if (line[0] == 'P')
    {
      nPlines++ ;
      sscanf(line, "P %s %s %s", auxstr[0], auxstr[1], auxstr[2]) ;
      for (k = 0 ; k < 3 ; k++)	gPBC[k] = gnumb(auxstr[k]) ;
    }
  }

//  // Check if everything was properly read:
//  fprintf(stderr, "%d %d %d\n", ngroups, nmols, natoms) ;
//  for (g = 1 ; g <= ngroups ; g++)
//  {
//    fprintf(stderr, "%s  %d\n", group[g].name, group[g].nmols) ;
//    for (i = 1 ; i <= group[g].nmols ; i++)
//    {
//      m = group[g].mol[i] ;
//      fprintf(stderr, "  %6d  %6d\n", mol[m].firstatom, mol[m].lastatom) ;
//    }
//  }
//  for (st = 1 ; st <= nstages ; st++)
//  {
//    fprintf(stderr, "O  %d\n", st) ;
//    for (o = nobjs[st-1]+1 ; o <= nobjs[st] ; o++)
//      fprintf(stderr, "O      %d\n", obj[o].mol) ;
//  }
//  fprintf(stderr, "C %d %d %d\n", gcent[0], gcent[1], gcent[2]) ;
//  fprintf(stderr, "P %d %d %d\n", gPBC[0], gPBC[1], gPBC[2]) ;

  // Some final checks for errors or warnings:
  for (g = 1 ; g <= ngroups ; g++)
    if (group[g].nmols == 0)
      message('E', "Group '%s' is empty.\n", group[g].name) ;
  for (a = 1 ; a <= natoms ; a++)
    if (atom[a].mol == 0)
      message('E', "Atom %d not assigned to a molecule.\n", a) ;
  if (nstages == 0) message('E', "Give at least one A-lines.\n") ;
  if (nClines != 1) message('E', "Give one (and only one) C-line.\n") ;
  if (nPlines != 1) message('E', "Give one (and only one) P-line.\n") ;

  fclose(fp_moldef) ;
}


void write_frame(void)
{
  int a ;

  printf("%s", title) ;  // already '\n'-terminated
  printf("%d\n", natoms) ;
  for (a = 1 ; a <= natoms ; a++)
  {
    printf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%s\n",
	   atom[a].rnumb, atom[a].rname, atom[a].aname, atom[a].anumb,
	   atom[a].r[0], atom[a].r[1], atom[a].r[2], atom[a].vels) ;
  }
  printf("%s", boxline) ;  // already '\n'-terminated
}


// Return the number of a group from its name.
int gnumb(char groupname[])
{
  int g, gg = -1 ;

  for (g = 0 ; g <= ngroups ; g++)
    if (strcmp(groupname, group[g].name) == 0) gg = g ;
  if (gg == -1) message('E', "Group '%s' not defined.\n", groupname) ;
  return gg ;
}


// Map from physical to scaled coordinates.
void r_to_s(void)
{
  int a, k, j ;

  for (a = 1 ; a <= natoms ; a++)
    for (k = 0 ; k < 3 ; k++)
    {
      atom[a].s[k] = 0 ;
      for (j = 0 ; j < 3 ; j++)
	atom[a].s[k] += iH[k][j]*(atom[a].r[j]-boxC[j]) ;
    }
}


// Map from scaled to physical coordinates.
void s_to_r(void)
{
  int a, k, j ;

  for (a = 1 ; a <= natoms ; a++)
    for (k = 0 ; k < 3 ; k++)
    {
      atom[a].r[k] = boxC[k] ;
      for (j = 0 ; j < 3 ; j++)
	atom[a].r[k] += H[k][j]*atom[a].s[j] ;
    }
}


// Compute inter-molecular distances as the one between their closest atoms.
void compute_objdist(void)
{
  int o1, o2, m1, m2, a1, a2, k, j ;
  float d2_min, drMI2, ds[3] = {0,0,0}, dsMI[3] = {0,0,0}, drMI[3] = {0,0,0} ;

  for (o1 = 1 ; o1 < nobjs_tot ; o1++)
  for (o2 = o1+1 ; o2 <= nobjs_tot ; o2++)
  {
    d2_min = BIG_DISTANCE ;
    m1 = obj[o1].mol ;
    m2 = obj[o2].mol ;
    for (a1 = mol[m1].firstatom ; a1 <= mol[m1].lastatom ; a1++)
    for (a2 = mol[m2].firstatom ; a2 <= mol[m2].lastatom ; a2++)
    {
      drMI2 = 0 ;
      for (k = 0 ; k < 3 ; k++)
      {
	ds[k] = atom[a2].s[k] - atom[a1].s[k] ;
	dsMI[k] = ds[k] - nint(ds[k]) ;
	drMI[k] = 0 ;
	for (j = 0 ; j < 3 ; j++) drMI[k] += H[k][j] * dsMI[j] ;
	drMI2 += drMI[k] * drMI[k] ;
      }
      if (drMI2 < d2_min)
      {
	dist[o1][o2].dRMI2 = dist[o2][o1].dRMI2 = d2_min = drMI2 ;
	for (k = 0 ; k < 3 ; k++)
	{
	  dist[o1][o2].dS[k] = ds[k] ;
	  dist[o2][o1].dS[k] = -ds[k] ;
	}
      }
    }
  }
}


// For N objects, this algorithm scales with N^3, but one scaling with N^2
// (e.g., based on clever distance sorting) can probably be devised.
void assemble_objs(void)
{
  int st, not_empty, o1, o2, o1_min = -1, o2_min = -1, a, k, m ;
  float d2_min ;

  for (o1 = 1 ; o1 <= nobjs_tot ; o1++) obj[o1].assembled = 0 ;

  // Use first object as nucleation point (doesn't matter which one is
  // used, since it always gives the same assembling, but translated):
  obj[1].assembled = 1 ;
  for (k = 0 ; k < 3 ; k++) obj[1].dN[k] = 0 ;

  // Cluster objects by sequential assembling stages:
  for (st = 1 ; st <= nstages ; st++)
  {
    do
    {
      d2_min = BIG_DISTANCE ;
      not_empty = 0 ;
      for (o1 = nobjs[st-1]+1 ; o1 <= nobjs[st] ; o1++)
      {
	if (obj[o1].assembled == 0)
	{
	  not_empty = 1 ;
	  for (o2 = 1 ; o2 <= nobjs[st] ; o2++)
	  {
	    // I tested switching the two conditions, but got no speed
	    // difference. But may be system-dependent...
	    if (obj[o2].assembled == 1)
	    {
	      if (dist[o1][o2].dRMI2 < d2_min)
	      {
		d2_min = dist[o1][o2].dRMI2 ;
		o1_min = o1 ;
		o2_min = o2 ;
	      }
	    }
	  }
	}
      }
      if (o1_min != -1 && o2_min != -1)
      {
	obj[o1_min].assembled = 1 ;
	for (k = 0 ; k < 3 ; k++)
	  // Since o2 is already assembled, its dN was already assigned.
	  obj[o1_min].dN[k] = obj[o2_min].dN[k] + nint(dist[o1_min][o2_min].dS[k]) ;
      }
    }
    while (not_empty) ;
  }

  // Apply box displacements from assembling:
  for (o1 = 1 ; o1 <= nobjs_tot ; o1++)
  {
    m = obj[o1].mol ;
    for (a = mol[m].firstatom ; a <= mol[m].lastatom ; a++)
      for (k = 0 ; k < 3 ; k++) atom[a].s[k] += obj[o1].dN[k] ;
  }
}


void center(void)
{
  int k, g, i, m, a ;
  float smin[3], smax[3], smid ;
  const char* ord[] = { "1st", "2nd", "3rd" } ;
  
  for (k = 0 ; k < 3 ; k++)
  {
    smin[k] = +BIG_DISTANCE ;
    smax[k] = -BIG_DISTANCE ;
    g = gcent[k] ;
    for (i = 1 ; i <= group[g].nmols ; i++)
    {
      m = group[g].mol[i] ;
      for (a = mol[m].firstatom ; a <= mol[m].lastatom ; a++)
      {
	if (atom[a].s[k] < smin[k]) smin[k] = atom[a].s[k] ;
	if (atom[a].s[k] > smax[k]) smax[k] = atom[a].s[k] ;
      }
    }
    smid = 0.5 * (smin[k] + smax[k]) ;
    for (a = 1 ; a <= natoms ; a++) atom[a].s[k] -= smid ;
    // Check if the box range was exceeded by the centering group and, if
    // it was, give an error or warning depending on centEW[k].
    if (smax[k] - smin[k] > 1)
    {
      message(centEW[k],
	      "Group '%s' exceeds range of %s box vector in frame %d.\n",
	      group[gcent[k]].name, ord[k], nframe) ;
    }
  }
}


void apply_PBCs(void)
{
  int k, g, i, m, a, dnCOM ;
  float sCOM ;

  for (k = 0 ; k < 3 ; k++)
  {
    g = gPBC[k] ;
    for (i = 1 ; i <= group[g].nmols ; i++)
    {
      m = group[g].mol[i] ;
      // Compute k-coordinate of m's "COM" (actually, geometric center):
      sCOM = 0 ;
      for (a = mol[m].firstatom ; a <= mol[m].lastatom ; a++)
	sCOM += atom[a].s[k] ;
      sCOM /= (mol[m].lastatom - mol[m].firstatom + 1) ;
      // Move m's atoms to inside the reference box:
      dnCOM = nint(sCOM) ;
      for (a = mol[m].firstatom ; a <= mol[m].lastatom ; a++)
	atom[a].s[k] -= dnCOM ;
    }
  }
}


// Fortran-like named. Useful property: nint(-x) = -nint(x)
int nint(float x)
{
  return floor(x + 0.5) ;
}


void message(char mtype, char *format, ...)
{
  va_list args ;
  const char cmd[] = "fixbox" ;
  const char usage[] = "Usage: %s GRO_FILE MOLDEF_FILE\nFixBox version 1.2\n" ;

  va_start(args, format) ;
  if (mtype != 'W' && mtype != 'E' && mtype != 'U')
    message('E', "Wrong use of message() function.\n") ;
  if (mtype == 'W') fprintf(stderr, "%s: WARNING: ", cmd) ;
  else fprintf(stderr, "%s: ERROR: ", cmd) ;
  vfprintf(stderr, format, args) ;
  va_end(args) ;
  if (mtype == 'U') fprintf(stderr, usage, cmd) ;
  if (mtype != 'W') exit(1) ;
}

