/*
 * main.cxx
 * 
 * Copyright 2024 Anirban <anirban.pal@protonmail.com>
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#include "headers.h"

int main(int argc, char **argv)
{
  //1. Initialize empty fiber network object and read data from file to create the network
  fibnetwork fn;
  fn.flog = fopen(LOG,"w");
  fprintf(fn.flog,"                 %6s %14s %14s %14s %14s %14s %14s %14s %14s\n","STEP", "AXIAL", "BENDING", "COHESIVE", "KINETIC", "TOTAL", "MOMENTUM_X", "MOMENTUM_Y", "MOMENTUM_Z");
	
  
  fn.readfile(argv[1]);
  
  //2. Print the initial structure of the fiber network, with coarse-grain control points and fine-grained form.
  fn.tstep = 0;
  fn.printlammps_cps((char*) DUMP1,(char*) "w");
  fn.printlammps((char*) DUMP2,(char*) "w", NBEZ);

  clock_t begin1,end1; 
  begin1 = clock(); 
  
  //3. Minimize the structure if needed. Perform numerical integration of equations of motion  
  loop(i,1) {
		tao_integrate(fn,NSTEP,dt);
		fn.minimize();
	}
  
  end1 = clock();
  double time1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("Total time taken(s) = %lf\n",time1);

	fclose(fn.flog);
  return 0;
}
