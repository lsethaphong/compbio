#!/usr/local/bin/perl

use Math::Trig;
####################################################
#
#  Parsing stdoutput files from the overlap
#  routine of Ptraj to get the cylindrical volumetric
#  density of ions within a z distance and radius r
#  of the precessing midpoint of a kissing loop
#  channel
#
#  9 June 2009, added bulk ion statistics
#               highest residency count
#  
#  30 May 2009, added selective surface plot
#  
#  19 May 2009, added surface plot feature
#     x,y,z -> z,r,ion_density 
#  May 2009, t_step is a variable
#
#  9 Sept 2009, changed interaction residence to 4 Angstroms
#               from 5 Angstroms
#
#  2 Oct 2009, sourced out ion pocket occupancy
#              to a data file *.pkt
#   
#  3 Jun 2010, added cumulative table and shortened
#              interaction distance from 4.0 
#
#  8 May 2012, added variable interaction distance
#              default is 5.0 Angstroms
#
#  Latsavongsakda Sethaphong
#  Yingling Lab 
#  March 2009
#
####################################################
# The cylinder construction equivalent
# to the construction of a parallelogram
#  @P1                              @P3
#   |              | -dz- |          |
#  X1--------------o---|_-----------X2
#   |               \     |          |
#  @P2               \    |         @P4
#                   r \   |(d_perp)
#                      \  |
#                       \ |
#                        \|
#                         X0
#
# 1. o, is the midpoint of the KL channel
# 2. X1, X2 are the midpoints between two phosphates
# 3. comprising the opening to the ion channel                        
# 4. d_perp is the perpendicular distance between point
# 5. X0 and the line intersecting X1 and X2
# 6. r is the ray to X0 from o.
# 7. dz = (||r||^2-||d_perp||^2)^(1/2)  
#                    
# Defining a line 
# given two points X1 and X2 
#
# a vector     |x1, + (x2-x1)t| 
# along the => |y1, + (y2-y1)t|
# line is      |z1, + (z2-z1)t|
#
# t_min is = -(X1-X0)dot(X2-X1)/||X2-X1||^2
#
# (AxB)^2 = A^2.B^2 - (A dot B)^2
#
# d_perp_squared = |X1-X0|^2.|X2-X1|^2-[(X1-X0)dot(X2-X1)]^2
#                  -----------------------------------------
#                                  |X2-X1|^2
# 
#
#                = |(X2-X1)x(X1-X0)|^2
#                  -------------------
#                        |X2-X1|^2
#
# d_perp         = |(X2-X1)x(X1-X0)|   ; one may reverse X0-X1 to reverse signage
#                  -----------------
#                       |X2-X1|
# Weisstein, Eric W. "Point-Line Distance--3-Dimensional."
#    From Mathworld -- a Wolfram Web Resource.
#    http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
#===================================================
# Defining the snaptshot PDB file format generated from
# ptraj
#
# ATOM ATOM_NUM TYPE RESNAME RESNUM X Y Z OCCUPANCY B-FACTOR
#  0    1        2   3         4    5 6 7     8         9 
# when performing a string split on spaces, one should expect these
# indexation numbers
#
#===================================================
# Input parameters then the atom numbers for the items of intererst
# <location of snapshots> 
# <phosphate 1 atom no> <phosphate 2 atom no> <phosphate 3> <phosphate 4>
# <atom number of first ion> <atom number of last ion>
# <output file name> 
# 
#===================================================
# Output format
# give the origin x y z
# give the R distance to origin
# give the d_z distance along the cylinder axis
# give the d_perp distance to the cylinder axis
# give the sodium ID

sub trim($);
# x,y,z operations
sub get_d_perp;
#sub get_r2x0;
sub get_z_dist;
sub get_midpt;
sub get_r_dist;
sub get_u_vec;
sub dot_prod;
sub cross_prod;
sub move_pt;
sub get_z_phi_idx;
sub get_z_phi_raw;
sub report_stats;

sub trim($)
{
   my $string=shift;
   $string =~ s/^\s+//;
   $string =~ s/\s+$//;
   return $string;
}

sub move_pt {
     my ($x0, $y0, $z0, $dx0, $dy0, $dz0, $mag_d) = @_;
     my @new_pt = (0.0,0.0,0.0);
#     print "$x0,$y0,$z0 moving to";
     $new_pt[0]= $x0 + $mag_d*$dx0;
     $new_pt[1]= $y0 + $mag_d*$dy0;
     $new_pt[2]= $z0 + $mag_d*$dz0;
#     print" @new_pt by $mag_d \n";
return @new_pt
}

#=====================================================
# Simple vector operations
# using the X1 as the common origin for both rays 
#
#-----------------------------------------------------
sub dot_prod {
    my ($x0,$y0,$z0,$x1,$y1,$z1,$x2,$y2,$z2) = @_; #take in the list of arguments
    my $d_perp = 0.0; # evaluated values
    my ($a1,$a2,$a3,$b1,$b2,$b3,$c1,$c2,$c3)= (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

# construct vector {X2-X1}
  $a1=($x2-$x1);#a1
  $a2=($y2-$y1);#a2
  $a3=($z2-$z1);#a3
#  $modulus_a =($a1**2+$a2**2+$a3**2)**0.5 ;
# construct vector {X1-X0}       
#  $b1=($x1-$x0);#b1
#  $b2=($y1-$y0);#b2
#  $b3=($z1-$z0);#b3
# take position 1 as the common origin
  $b1=($x0-$x1);#b1
  $b2=($y0-$y1);#b2
  $b3=($z0-$z1);#b3

#  $modulus_b =($b1**2+$b2**2+$b3**2)**0.5 ;
# construct the dot product of {X2-X1}x{X1-X0}
  $c1=$a1*$b1;
  $c2=$a2*$b2;
  $c3=$a3*$b3;
#  $modulus_c = ($c1**2+$c2**2+$c3**2)**0.5;
  $dot_value = $c1+$c2+$c3;
return $dot_value
}
#=====================================================
# Generate the
#
#
#=====================================================
# This function returns the perpendicular distance from
# a point to a line defined by X1 and X2 in 3-D space
# Call the subroutine with an & e.g. $val=&myroutine
# The data will be read as strings from the file, and
# must be converted into a floating point/double number
#-----------------------------------------------------
sub cross_prod  {
#=====================================================
# X0 X1 and X2 are inputs, to generate the length d_perp
# c.f. construction shown in the beginning
#-----------------------------------------------------
   my ($x0,$y0,$z0,$x1,$y1,$z1,$x2,$y2,$z2) = @_; #take in the list of arguments
   my $d_perp = 0.0; # evaluated values
   my ($a1,$a2,$a3,$b1,$b2,$b3,$c1,$c2,$c3)= (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
   my ($c1_n, $c2_n, $c3_n) = (0.0,0.0,0.0);
   my @xx_vec = (0.0,0.0,0.0,0.0,0.0,0.0,0.0); # cross vector product with normalized unit  
# construct vector {X2-X1}
   $a1=($x2-$x1);#a1
   $a2=($y2-$y1);#a2
   $a3=($z2-$z1);#a3
   $modulus_a =($a1**2+$a2**2+$a3**2)**0.5;
# construct vector {X1-X0}       
   $b1= ($x0-$x1);#b1
   $b2= ($y0-$y1);#b2
   $b3= ($z0-$z1);#b3
   $modulus_b =($b1**2+$b2**2+$b3**2)**0.5;
# construct the cross product of {X2-X1}x{X1-X0}
   $c1=$a2*$b3-$a3*$b2;
   $c2=$a3*$b1-$a1*$b3;
   $c3=$a1*$b2-$a2*$b1;
   $modulus_c = ($c1**2+$c2**2+$c3**2)**0.5; # area of parallelogram formed by rays X10 and X1
   $d_perp = $modulus_c/$modulus_a;
#normalized unit vector
   $c1_n=$c1/$modulus_c; 
   $c2_n=$c2/$modulus_c;
   $c3_n=$c3/$modulus_c;
   @xx_vec = ($c1, $c2, $c3, $d_perp, $c1_n, $c2_n, $c3_n);
return @xx_vec;
}

#====================================================
# Calculates the magnitude of a ray, r, pointing
# from position X1 to position X2
#----------------------------------------------------
sub get_r_dist {
    my ($x1, $y1, $z1, $x2, $y2, $z2) = @_;
    my $r_dist = 0;
    my ($a1,$a2,$a3) = (0.0,0.0,0.0);
    my $modulus_a = 0.0;
# construct vector {X2-X1}
  $a1=($x2-$x1);#a1
  $a2=($y2-$y1);#a2
  $a3=($z2-$z1);#a3
  $modulus_a =($a1**2 +$a2**2 + $a3**2)**0.5 ;

# return the x,y,z and modulus
  @r_dist= ($a1,$a2,$a3,$modulus_a);
return @r_dist;
}

#======================================================
# This routine returns the mid point between X1 and X2
# position X1 points to position X2
# the magnitude of this ray is also computed
#------------------------------------------------------
sub get_midpt {
   my ($x1,$y1,$z1,$x2,$y2,$z2) = @_;
   my ($xm,$ym,$zm) = (0.0,0.0,0.0);
   my ($xp,$yp,$zp) = (0.0,0.0,0.0);
   my @mid_point = (0.0,0.0,0.0); # only 
   my $rm = 0.0;

   $xm= ($x2-$x1);
   $ym= ($y2-$y1);
   $zm= ($z2-$z1);
   $rm= ($xm**2 + $ym**2 + $zm**2)**0.5; # radial distance

   $xp= 0.5*$xm + $x1;
   $yp= 0.5*$ym + $y1;
   $zp= 0.5*$zm + $z1;
   @mid_point = ($xp,$yp,$zp);#,$rm);
return  @mid_point; # returns the x,y,z,modulus of the midpoint
}
#====================================================
# This routine normalizes the ray pointing from x1 to x2
# as a unit vector. The magnitude of the
# actualray is given as the third component
#-----------------------------------------------------
sub get_u_vec {
   my ($x1,$y1,$z1,$x2,$y2,$z2) = @_;
   my ($xr,$yr,$zr) = (0.0,0.0,0.0);
   my ($xp,$yp,$zp) = (0.0,0.0,0.0);
   my @u_vec= (0.0,0.0,0.0,0.0);
   $xr= ($x2-$x1);
   $yr= ($y2-$y1);
   $zr= ($z2-$z1);
   $rr= ($xr**2 + $yr**2 + $zr**2)**0.5; # radial distance
   
#normalize the unit vector
   $xr=$xr/$rr;
   $yr=$yr/$rr;
   $zr=$zr/$rr;
#   $xp= 0.5*$xm + $x1;
#   $yp= 0.5*$ym + $y1;
#   $zp= 0.5*$zm + $z1;
   @u_vec = ($xr,$yr,$zr,$rr);
return @u_vec;
}

#======================================================
# Returns the geometric difference between two numbers 
#
#------------------------------------------------------
sub get_z_dist {
   
   my ($c_z,$b_z) = @_;
   my $d_z = 0.0; # float
   if ($c_z > $b_z) {
      $d_z=($c_z**2-$b_z**2)**0.5;
   } else {
      $d_z=($b_z**2-$c_z**2)**0.5;
   }
   
return $d_z;
}

#======================================================
# Returns the z index of a given r,z pair value
# r,z are incremented in 0.5 blocks
# z is in +/- the input value
# may make the block size arbitrary later
#------------------------------------------------------
sub get_zr_index {
   my ($z_raw, $r_raw, $c_val) = @_;
   my $zr_idx = 0; # integer 
   my ($z_idx, $r_idx) = (0, 0);
              
   if ($c_val > 0) {       # $z_max is defined globally
      $z_idx = int (($z_raw + $z_max)/0.5); # x form value
   } else {
      $z_idx = int (($z_max - $z_raw)/0.5); # x form value
   }

   $r_idx = int ($r_raw/0.5); # y form value

   $zr_idx = ($z_idx)%(($z_max/0.5)*2) + ($r_idx*(($z_max/0.5)*2));
return $zr_idx;
}
#======================================================
# Returns the raw z and r grid values in order to
# generate the x,y,z surface map coordinates
# row major order
#------------------------------------------------------
sub get_zr_raw {

    my ($zr_idx) = @_;
    my @zr_raw = (0.0,0.0);
    my ($z_idx,$r_idx)= (0.0,0.0);

#   print "raw index is \ $zr_idx \n";    

   $z_idx = ($zr_idx)%(($z_max/(0.5))*2);
   $r_idx = $zr_idx - $z_idx;
   $r_idx = (($r_idx/(($z_max/(0.5))*2))*0.5); # converting back to half units

   $z_idx = ((($z_idx)-(2*$z_max))*0.5);
   @zr_raw = ($z_idx, $r_idx);

    
return @zr_raw;
}

#======================================================
# Z and Phi tabulation routines
#------------------------------------------------------
# phi is fed in radians, but will be stored as degrees
#------------------------------------------------------
sub get_z_phi_idx {
   my ($z_raw, $phi_raw,) = @_;
   my $z_phi_idx = 0; # integer 
   my ($z_idx, $phi_idx) = (0, 0);
   # z_raw is fed signed
   $z_idx = int (($z_raw + $z_max)/0.5); # x form value
   #phi_raw is fed signed 
   $pi = 3.1415926535897932384626433832795;
   if ($phi_raw > 0) {
     $phi_idx = int ( ((($phi_raw/$pi)*180))/5); # y form value 0 to 76
   } else {
     $phi_idx = int ( ((360 + ($phi_raw/$pi)*180) )/5); # y form value 0 to 76
   }
   print "phi_idx $phi_idx \n";
   $z_phi_idx = ($z_idx)%(($z_max/0.5)*2) + ($phi_idx*(($z_max/0.5)*2));


return $z_phi_idx;
}

sub get_z_phi_raw {

    my ($z_phi_idx) = @_;
    my @z_phi_raw = (0.0,0.0);
    my ($z_idx,$phi_idx)= (0.0,0.0);

#   print "raw index is \ $zr_idx \n";    

   $z_idx = ($z_phi_idx)%(($z_max/(0.5))*2); # module remainder
   $phi_idx = $z_phi_idx - $z_idx;#0 to 71 # remainder in x/phi direction

#   if ($phi_idx < 36) {
#      $phi_idx = ($phi_idx*5); # converting back to half units
#   } elsif ($phi_idx > 36) { 
#      $phi_idx = $phi_idx*5 - 360; 
#   } else {
#      $phi_idx = 180; # equivalent to zero
#   }
      $phi_idx = ($phi_idx*5)/(($z_max/0.5)*2);

   $z_idx = ((($z_idx)-(2*$z_max))*0.5);

   @z_phi_raw = ($z_idx, $phi_idx);

return @z_phi_raw;
}

# write statistics and data files
sub report_stats {

                 open(MYFILE, '>>'.$base.'.'.'out')|| die;
        #=======================================
        # calculate occupancy statistics
        #---------------------------------------
                 $mean_occu = 0; 
                 $std_dev_occu = 0;
                 $mean_occu = ($ion_total/$frame_cnt);
                 $variance_occu = ($sqr_occu_sum/$frame_cnt) - $mean_occu*$mean_occu;
                 $std_dev_occu = sqrt($variance_occu);

         print MYFILE "Occupancy Statistics -- var= \ $variance_occu \ \ mean= \ $mean_occu \ \ +/- \ \ $std_dev_occu  \n";

         #        close (MYFILE);

         print "Occupancy Statistics -- var= \ $variance_occu \ \ mean= \ $mean_occu \ \ +/- \ \ $std_dev_occu  \n";
         print "Totals= \ \ $ion_total  \ \n";
         #======================================
         # calculate bulk statistics 
         #--------------------------------------
                 $mean_bulk = 0;
                 $std_dev_bulk = 0;
                 $mean_bulk = ($bulk_total/$frame_cnt);
                 $variance_bulk = ($sqr_bulk_sum/$frame_cnt) - $mean_bulk*$mean_bulk;
                 $std_dev_bulk = sqrt($variance_bulk);

         print MYFILE "Bulk Occupancy Statistics -- var= \ $variance_bulk \ \ mean= \ $mean_bulk \ \ +/- \ \ $std_dev_bulk \n";

        #         close (MYFILE);
         print "Bulk Occupancy Statistics -- var= \ $variance_bulk \ \ mean= \ $mean_bulk \ \ +/- \ \ $std_dev_bulk \n";
         print "Totals= \ \ $bulk_total \ \n";
                 $i = 0;
                 $i_dum = 0;
                 $highest_time = 0;
                 $highest_id = 0;
                 while ($i < $ion_high_numid) {
                  $i_dum = $i + 1;
                  if ($ion_high_residency_cnt[$i_dum] > $highest_time) {
                     $highest_time = $ion_high_residency_cnt[$i_dum];
                     $highest_id = $i_dum;
                  }

                 # print "Ion -- $i_dum \ \ No.   \ \ $ion_residency_cnt[$i_dum] \ \n";
                  $i = $i + 1;
                  }
                 $highest_time = $t_step*$highest_time;
                 print "Highest Occupancy Ion -- \ $highest_id \ \ time= \ $highest_time \ \n";

         #======================================
         # calculate cylinder statistics
         #--------------------------------------
                 $mean_cylinder = 0;
                 $std_dev_cylinder = 0;
                 $mean_cylinder = ($cylinder_total/$frame_cnt);
                 $variance_cylinder = ($sqr_cylinder_sum/$frame_cnt) - $mean_cylinder*$mean_cylinder;
                 $std_dev_cylinder = sqrt($variance_cylinder);

         print MYFILE "Cylinder Occupancy Statistics -- var= \ $variance_cylinder \ \ mean= \ $mean_cylinder \ \ +/- \ \ $std_dev_cylinder \n";

                 close (MYFILE);

         print "Cylinder Occupancy Statistics -- var= \ $variance_cylinder \ \ mean= \ $mean_cylinder \ \ +/- \ \ $std_dev_cylinder \n";
         print "Totals= \ \ $cylinder_total \ \n";
                 $i = 0;
                 $i_dum = 0;
                 $highest_time = 0;
                 $highest_id = 0;
                 while ($i < $ion_high_numid) {
                  $i_dum = $i + 1;
                  if ($ion_cylinder_high_residency_cnt[$i_dum] > $highest_time) {
                     $highest_time = $ion_cylinder_high_residency_cnt[$i_dum];
                     $highest_id = $i_dum;
                  } 

                 # print "Ion -- $i_dum \ \ No.   \ \ $ion_cylinder_residency_cnt[$i_dum] \ \n";
                  $i = $i + 1; 
                  }
                 $highest_time = $t_step*$highest_time;
                 print "Highest Cylinder Occupancy Ion -- \ $highest_id \ \ time= \ $highest_time \ \n";

         #======================================
         # dump ion grid count
         #--------------------------------------
                 open(MYFILE, '>'.$base.'.'.'grid');
                 close(MYFILE);

                 open(MYFILE, '>>'.$base.'.'.'grid')||die;
                 print "grid size is \ $ion_grid_sz \n";

                 $g_idx = 0;

                 for ($grid_idx = 0; $grid_idx < $ion_grid_sz; $grid_idx++) {
                    $g_idx = $grid_idx + 1; # taking advantage
                    @zr_raw = get_zr_raw($g_idx);
                 #   print "$zr_raw[0] \ $zr_raw[1] \ $ion_grid_cnt[$g_idx] \n";
                    if ($ion_grid_cnt[$g_idx] > 0) {
                       print MYFILE "$zr_raw[0] \ $zr_raw[1] \ $ion_grid_cnt[$g_idx] \ \n";
                    }
                 }

                 close(MYFILE);

            #===================================
            # dump residue HIGHEST residence ion times into tabular format
            #-----------------------------------
                 open(MYFILE, '>'.$base.'.'.'hire');
                 close(MYFILE);

                 open(MYFILE, '>>'.$base.'.'.'hire')||die;
                print "writing ion high occupancies \n";

                 $g_idx = 0;

                 print MYFILE "ION ByRes ";
                 for ($ii_idx = 0; $ii_idx < $res_B; $ii_idx++) {
                       $iii_idx = $ii_idx + 1;
                       print MYFILE ", $iii_idx ";
                       }
                 print MYFILE "\n";
                 
                 for ($ion_idx = 0; $ion_idx < $bulk_ion_tot; $ion_idx++) {
                       $is_ion = $ion_idx + 1;
#                      print "doing $res_B residues for ion $is_ion\n";
                      for ($res_idx = 0; $res_idx < $res_B; $res_idx++)  {
                         $is_res = $res_idx + 1;
                         
                         $g_idx = ($ion_idx )*$res_B + $is_res;

                         if ($res_idx == 0) {
 #                           print "dump one  at $g_idx \n";
                            if ($ion_is_here_high[$g_idx] > 0) {
                              print MYFILE " $is_ion , \ $ion_is_here_high[$g_idx] \ ";
                            } else {
                              print MYFILE " $is_ion, 0 \ ";
                            }
                         } else {
                            if ($ion_is_here_high[$g_idx] > 0) {
                              print MYFILE " , \ $ion_is_here_high[$g_idx] \ ";
                            } else {
                              print MYFILE " , 0 \ ";
                            }
                         }
                       }
                   
                      print MYFILE "\n"; # end of line
  #                    print "finished ion $is_ion \n";
                 }

                 close(MYFILE);

             #===================================
             # HIGHEST Ion Residence Times Dump in x,y,z format for tecplot
             #-----------------------------------
                 open(MYFILE, '>'.$base.'.'.'hitc');
                 close(MYFILE);

                 open(MYFILE, '>>'.$base.'.'.'hitc')||die;
                print "writing ion high occupancies for tecplot \n";

                 $g_idx = 0;

                 for ($ion_idx = 0; $ion_idx < $bulk_ion_tot; $ion_idx++) {
                       $is_ion = $ion_idx + 1;
#                      print "doing $res_B residues for ion $is_ion\n";
                      for ($res_idx = 0; $res_idx < $res_B; $res_idx++)  {
                         $is_res = $res_idx + 1;

                         $g_idx = ($ion_idx )*$res_B + $is_res;

 #                           dumping by ion as x, y as residue, occupancy as z
                            if ($ion_is_here_high[$g_idx] > 0) {
                              print MYFILE " $is_ion, $is_res, $ion_is_here_high[$g_idx] \n";
                            }# else {
                             # print MYFILE " $is_ion, $is_res, 0 \n";
                           # }
                          
                       }

                 }

                 close(MYFILE);
             #-----------------------------------
             #===========================================================
             # DUMPING CUMULATIVE DATA
             #===========================================================
            #===================================
            # CUMULATIVE Residence ion times into tabular format
            #-----------------------------------
                 open(MYFILE, '>'.$base.'.'.'tore');
                 close(MYFILE);
                 open(MYFILE, '>>'.$base.'.'.'tore')||die;
                 print "writing cumulative ion occupancies \n";

                 $g_idx = 0;

                 print MYFILE "ION ByRes ";
                 for ($ii_idx = 0; $ii_idx < $res_B; $ii_idx++) {
                       $iii_idx = $ii_idx + 1;
                       print MYFILE ", $iii_idx ";
                       }
                 print MYFILE "\n";

                 for ($ion_idx = 0; $ion_idx < $bulk_ion_tot; $ion_idx++) {
                       $is_ion = $ion_idx + 1;
#                      print "doing $res_B residues for ion $is_ion\n";
                      for ($res_idx = 0; $res_idx < $res_B; $res_idx++)  {
                         $is_res = $res_idx + 1;

                         $g_idx = ($ion_idx )*$res_B + $is_res;

                         if ($res_idx == 0) {
 #                           print "dump one  at $g_idx \n";
                            if ($ion_is_here_tot[$g_idx] > 0) {
                              print MYFILE " $is_ion , \ $ion_is_here_tot[$g_idx] \ ";
                            } else {
                              print MYFILE " $is_ion, 0 \ ";
                            }
                         } else {
                            if ($ion_is_here_high[$g_idx] > 0) {
                              print MYFILE " , \ $ion_is_here_tot[$g_idx] \ ";
                            } else {
                              print MYFILE " , 0 \ ";
                            }
                         }
                       }

                      print MYFILE "\n"; # end of line
  #                    print "finished ion $is_ion \n";
                 }

                 close(MYFILE);
             #------------------------------------

             #===================================
             # HIGHEST Ion Residence Times Dump in x,y,z format for tecplot
             #-----------------------------------
                 open(MYFILE, '>'.$base.'.'.'totc');
                 close(MYFILE);

                 open(MYFILE, '>>'.$base.'.'.'totc')||die;
                print "writing cumulative ion occupancies for tecplot \n";

                 $g_idx = 0;

                 for ($ion_idx = 0; $ion_idx < $bulk_ion_tot; $ion_idx++) {
                       $is_ion = $ion_idx + 1;
#                      print "doing $res_B residues for ion $is_ion\n";
                      for ($res_idx = 0; $res_idx < $res_B; $res_idx++)  {
                         $is_res = $res_idx + 1;

                         $g_idx = ($ion_idx )*$res_B + $is_res;

 #                           dumping by ion as x, y as residue, occupancy as z
                            if ($ion_is_here_tot[$g_idx] > 0) {
                              print MYFILE " $is_ion, $is_res, $ion_is_here_tot[$g_idx] \n";
                            }# else {
                             # print MYFILE " $is_ion, $is_res, 0 \n";
                           # }

                       }

                 }

                 close(MYFILE);
             #-----------------------------------

             #===================================
             # cumulative frequency by residue
             # dump ion tabular format
             #===================================
                 open(MYFILE, '>'.$base.'.'.'tab');
                 close(MYFILE);

                 open(MYFILE, '>>'.$base.'.'.'tab')||die;
                 $ion_grid_sz = $z_max*4;
                 print "grid size is \ $ion_grid_sz \n";

                 $g_idx = 0;
                    
                    print MYFILE "z/r"; # first entry 
                   # generate header
                   for ($grid_idx = 0; $grid_idx < $z_max*4; $grid_idx++) {

                      @zr_raw = get_zr_raw($grid_idx);
                      print MYFILE "\ @zr_raw[0]";
                    }
                    print MYFILE " \n";

                 $g_idx = 0;
                 for ($grid_idx = 0; $grid_idx < $ion_grid_sz; $grid_idx++) {
                    $g_idx = $grid_idx + 1; # taking advantage
                    @zr_raw = get_zr_raw($g_idx);

#                    print "$zr_raw[0] \ $zr_raw[1] \ $ion_grid_cnt[$g_idx]";
                    if ($ion_grid_cnt[$g_idx] > 0) {

                      if ($zr_raw[0] < -5.5) {
                        print MYFILE "$zr_raw[1] \ $ion_grid_cnt[$g_idx]";
                      } else {
                        if ($zr_raw[0] < 5.5) {
                          print MYFILE "\ $ion_grid_cnt[$g_idx]";
                        } else {
                          print MYFILE "\ $ion_grid_cnt[$g_idx] \n";
                        }
                      }

                    } else {

                      if ($zr_raw[0] < -5.5) {
                        print MYFILE "$zr_raw[1] \ 0";
                      } else {
                        if ($zr_raw[0] < 5.5) {
                          print MYFILE "\ 0";
                        } else {
                          print MYFILE "\ 0 \n";
                        }
                      }

                    }

                 }

                 close(MYFILE);

         #======================================
         # dump z phi count
         #--------------------------------------

         # dump residue count
                 open(MYFILE, '>'.$base.'.'.'zphi');
                 close(MYFILE);

                 open(MYFILE, '>>'.$base.'.'.'zphi')||die;
                 #print "grid size is \ $ion_grid_sz \n";

                 $g_idx = 0;
                 $z_phi_grid_sz = ($z_max/0.5)*2*(360/5);
                 print "grid size is $z_phi_grid_sz \n";
                 for ($grid_idx = 0; $grid_idx < $z_phi_grid_sz; $grid_idx++) {
                    $g_idx = $grid_idx + 1; # taking advantage
                 #   $idx_z_phi = get_z_phi_idx($dzc, $phi_value); # increment phi by 0.0872664625997 radians or 5 degrees
                 #   $z_phi_table[$idx_z_phi]=$z_phi_table[$idx_z_phi] + 1;
                    #print "made z phi index, $idx_z_phi for $dzc, $phi_value\n";
                    @z_phi_raw = get_z_phi_raw($g_idx);
                    if ($z_phi_table[$g_idx] > $thresh ) { # was 0
                 #   print "$zr_raw[0] \ $zr_raw[1] \ $ion_grid_cnt[$g_idx] \n";
                    print MYFILE "$z_phi_raw[0] \ $z_phi_raw[1] \ $z_phi_table[$g_idx] \ \n";
                    }
                 }

                 close(MYFILE);


         # Dump stemp B z phi table
                 open(MYFILE, '>'.$base.'.'.'zphiA');
                 close(MYFILE);

                 open(MYFILE, '>>'.$base.'.'.'zphiA')||die;
                 #print "grid size is \ $ion_grid_sz \n";

                 $g_idx = 0;
                 $z_phi_grid_sz = ($z_max/0.5)*2*(360/5);
                 print "grid sizeA is $z_phi_grid_sz \n";
                 for ($grid_idx = 0; $grid_idx < $z_phi_grid_sz; $grid_idx++) {
                    $g_idx = $grid_idx + 1; # taking advantage
                 #   $idx_z_phi = get_z_phi_idx($dzc, $phi_value); # increment phi by 0.0872664625997 radians or 5 degrees
                 #   $z_phi_table[$idx_z_phi]=$z_phi_table[$idx_z_phi] + 1;
                    #print "made z phi index, $idx_z_phi for $dzc, $phi_value\n";
                    @z_phi_raw = get_z_phi_raw($g_idx);
                    if ($z_phi_tableA[$g_idx] > $thresh) { # was zero
                 #   print "$zr_raw[0] \ $zr_raw[1] \ $ion_grid_cnt[$g_idx] \n";
                    print MYFILE "$z_phi_raw[0] \ $z_phi_raw[1] \ $z_phi_tableA[$g_idx] \ \n";
                    }
                 }

                 close(MYFILE);

         # Dump stem B z phi table
                 open(MYFILE, '>'.$base.'.'.'zphiB');
                 close(MYFILE);

                 open(MYFILE, '>>'.$base.'.'.'zphiB')||die;
                 #print "grid size is \ $ion_grid_sz \n";

                 $g_idx = 0;
                 $z_phi_grid_sz = ($z_max/0.5)*2*(360/5);
                 print "grid sizeB is $z_phi_grid_sz \n";
                 for ($grid_idx = 0; $grid_idx < $z_phi_grid_sz; $grid_idx++) {
                    $g_idx = $grid_idx + 1; # taking advantage
                 #   $idx_z_phi = get_z_phi_idx($dzc, $phi_value); # increment phi by 0.0872664625997 radians or 5 degrees
                 #   $z_phi_table[$idx_z_phi]=$z_phi_table[$idx_z_phi] + 1;
                    #print "made z phi index, $idx_z_phi for $dzc, $phi_value\n";
                    @z_phi_raw = get_z_phi_raw($g_idx);
                    if ($z_phi_tableB[$g_idx] > $thresh ) { # was zero
                 #   print "$zr_raw[0] \ $zr_raw[1] \ $ion_grid_cnt[$g_idx] \n";
                    print MYFILE "$z_phi_raw[0] \ $z_phi_raw[1] \ $z_phi_tableB[$g_idx] \ \n";
                    }
                 }

                 close(MYFILE);


                 
                 #===================================
                 # Base Frequency tables and close neighbor statistics
                 #-----------------------------------
                 open(MYFILE, '>'.$base.'.'.'resf');
                 close(MYFILE);

                 $g_idx = 0;
                 $my_type = -1;
                 $t_idx = 0;
                 open(MYFILE, '>>'.$base.'.'.'resf')||die;
                 for ($grid_idx = 0; $grid_idx < $res_max ; $grid_idx++)    {
                      $g_idx = $grid_idx + 1;
                      $t_idx = 7*$grid_idx; # 0 to 6 for seven cells
                           
                      $my_type = $atom_type_stat[$t_idx + 6]; # get residue type
                      
                      if ($residue_box[$g_idx] > 0) {
                         print MYFILE "$g_idx $residue_box[$g_idx] ";}
                      else {
                         print MYFILE "$g_idx 0 ";}
                      print "$my_type";
                      if    ($my_type > 0 && $my_type < 2) { print MYFILE " G ";}
                      elsif ($my_type > 1 && $my_type < 3) { print MYFILE " A ";}
                      elsif ($my_type > 2 && $my_type < 4) { print MYFILE " C ";}
                      elsif ($my_type > 3 && $my_type < 5) { print MYFILE " U ";}
                      elsif ($my_type =~ m/0/) { print MYFILE " T ";} else {print MYFILE " | "; print " | \n"; }
                      
                      $total_type = $atom_type_stat[$t_idx + 2];
                      if ($total_type > 0) {
                         $mean_type =   $atom_type_stat[$t_idx + 1]/$total_type;
                         $sqr_type_sum = $atom_type_stat[$t_idx + 0];
                      
                         $std_dev_type = 0;
                        
                         $variance_type = ($sqr_type_sum/$total_type) - $mean_type*$mean_type;
                         $std_dev_type = sqrt($variance_type);

                         if ($my_type > 0 && $my_type < 2) {    # G
                            print MYFILE " $total_type O6 \ $mean_type +/- $std_dev_type ";
                         } elsif ($my_type < 1) { # T
                            print MYFILE " $total_type O4 \ $mean_type +/- $std_dev_type ";
                         } elsif ($my_type > 2) { # C U
                            print MYFILE " $total_type O4 \ $mean_type +/- $std_dev_type ";
                         } else { 
                            print MYFILE " 0 NA \ 0  +/- 0 "; # not applicable for A
                         }

                       }

                      $total_type = 0;
                      $total_type = $atom_type_stat[$t_idx + 5];
                      if ($total_type > 0) {
                         $mean_type =   $atom_type_stat[$t_idx + 4]/$total_type;
                         $sqr_type_sum = $atom_type_stat[$t_idx + 3];
                      
                         $std_dev_type = 0;
                        
                         $variance_type = ($sqr_type_sum/$total_type) - $mean_type*$mean_type;
                         $std_dev_type = sqrt($variance_type);

                         if($my_type > 0 && $my_type < 2) { #G
                            print MYFILE " $total_type N7 $mean_type +/- $std_dev_type ";
                          } elsif($my_type > 1 && $my_type < 3) { #A
                            print MYFILE " $total_type N7 $mean_type +/- $std_dev_type ";
                          } else { #T,C,U
                            print MYFILE " $total_type O2 $mean_type +/- $std_dev_type ";
                          }
                       }


 
                      print MYFILE "\n";

                   }
                 close(MYFILE);

                 print " @residue_type \n";

return 0;
}
#======================================================
# END OF SUBROUTINE DECLARATIONS
#======================================================
#  BEGIN MAIN PROGRAM 
#======================================================
# if the last element is les than zero, no arguments were
# passed
if($#ARGV < 0){
   print "  Incorrect usage....\n";
   exit;
}


   $j=1;
   $atom_cnt=0; # number of atoms per set in proximity
   $frame_cnt=0; # frame occurrence
   $atom_total=0; # number of all atoms in the area for all frames
   $k=1;

   $base =     $ARGV[0]; # base name of pdb files
   $phos1=     $ARGV[1]*1.0; # Atom number of the first phosphorus
   $phos2=     $ARGV[2]*1.0; # Atom number of the second phosphorus
   $phos3=     $ARGV[3]*1.0; # Atom number of the third phosphorus
   $phos4=     $ARGV[4]*1.0; # Atom number of the fourth phosphorus
   $ion_start= $ARGV[5]*1.0; # First Ion's Atom Number
   $ion_stop=  $ARGV[6]*1.0; # Last Ion Atom Number
   $r_max=     $ARGV[7]*1.0; # Max radius of cylinder
   $z_max=     $ARGV[8]*1.0; # Max Length of cylinder from center
   $t_step=    $ARGV[9]*1.0; # expected time between frames
   $t_max =    $ARGV[10]*1.0; # optional max length of time
   $br_min =   $ARGV[11]*1.0; # radius of RNA plus 5 Angstroms
   $bz_min =   $ARGV[12]*1.0; # 1/2 z length of entire RNA
   $pi_n   =   $ARGV[13]*1.0; # pi division by n used to check helical ion distribution
   $t_start =  $ARGV[14]*1.0; # starting time
   $res_A =    $ARGV[15]*1.0; # end of strand A
   $res_B =    $ARGV[16]*1.0; # end of strand B 
   $thresh =   $ARGV[17]*1.0; # threshold occupancy for z phi charts
   $ia_dist =  $ARGV[18]*1.0; # minimum interactin distance btw ion & nearest nuclei
#   $res_T =    $ARGV[18]*1.0; # occupancy of specific residue id

  print  "p1\ $phos1 \ p2\  $phos2 \ p3\  $phos3 \ p4\ $phos4\n";
  print  "Na start \ $ion_start \ Na stop \ $ion_stop \n";
  print  "Channel radial dist max \ $r_max \ radial axis max \ $z_max \n";
  print  "Bulk radial dist min \ $br_min \ axis min \ $bz_min \n";
  print  "Residues Number \ $res_B \ min occupancy is \ $thresh \n"; 
  print  "start time is \ $t_start \ stop time is \ $t_max \n";
#foreach $i ( 0.. $#ARGV ) {
   open(MYFILE, '>'.$base.'.'.'out')|| die; # overwrite old files
   close(MYFILE);
   $k = 1;
   if ($t_start > 0) {
   $k =$t_start;
   }

# atomic interaction distance
  if ($ia_dist < 1.0) {
    $ia_dist = 3.2; # default interaction distance
   }

  print "Interaction distance is \ $ia_dist \n";

   if ($pi_n > 0) {
   $pi_n = 180/$pi_n; # in degrees 
   }
   $frame_cnt = 0;
   $ion_total = 0;
   $atom_numid = 0;
   $ion_numid = 0;
   $attempfile = 0;
# do for the first 20 ns
  $bulk_ion_tot = $ion_stop - $ion_start + 1.0; # total bulk ions

  print "total ions \ $bulk_ion_tot \ \n";
  print  "--------- Beginning Analysis ----------------------- \n";

# allow for proper counting
   $ion_start = $ion_start - 1;
   $ion_stop = $ion_stop + 1;

#
    
    $r_time = 400;
  if ($t_step > 1) {
    $r_time = $t_step + 20000/$t_step;
    if ($t_max > 1) {
       $r_time = $t_step + ($t_max - $k)/$t_step;
       #print "$r_time\n";
    }
  } 
#============================================================
# Statistical Variables 
#============================================================
# Circular distribution by phi
#------------------------------------------------------------
  @atom_tab=(0,0,0,0,0); #initial $d_idx
  $d_idx = 0;
  @residue_box=(0,0,0); # will grow to res_max
  @residue_type=(0,0,0,0,0); # only use the last four
#============================================================ 
# channel occupancy
#------------------------------------------------------------
  $sqr_occu_sum = 0;
  $occu_cnt = 0;
  $variance_occu =0;

#============================================================
# bulk occupancy statistics
#------------------------------------------------------------
#  $bulk_ion_tot = $ion_stop - $ion_stop + 1.0; # total bulk ions:w
  $sqr_bulk_sum  = 0;
  $bulk_cnt   = 0; # observed ions in the bulk per frame
  $bulk_total = 0; # total ions observed in bulk
  $variance_bulk = 0;

  $ion_high_numid = 0;
  @ion_residency_cnt = (0);
  @ion_high_residency_cnt = (0);
  
  $ion_num_idx = 0; # index of ion for the following arrays
  @ion_is_here_flag = (0); # flag for ion being present
  @ion_is_here_cnt = (0);
  @ion_is_here_high = (0); # residency near a residue
  @ion_is_here_tot = (0); # total cumulative residency
  $ion_is_here_dist = 3.2; # distance an ion is considered being present
#============================================================
# cylinder occupancy statistics
#------------------------------------------------------------
  $sqr_cylinder_sum  = 0;
  $cylinder_cnt   = 0; # observed ions in the bulk per frame
  $cylinder_total = 0; # total ions observed in bulk
  $variance_cylinder = 0;

#  $ion_high_numid = 0;
  @ion_cylinder_residency_cnt = (0);
  @ion_cylinder_high_residency_cnt = (0);


#============================================================
# surface topology array counter
#------------------------------------------------------------
   $zr_index = 0;
   @ion_grid_cnt = (0.0);
   $ion_grid_sz = 0; # initialize the grid size values
   $ion_grid_sz = int (($z_max*4)*($r_max/0.5));
   $ion_grid_cnt[$ion_grid_sz] = 0.0; # size it up to 

#=============================================================
# Residue by Atom Array Bins
#-------------------------------------------------------------
  @atom_type_name=(0); # G:1,2-> O6,N7; A:3->N7
  @atom_type_stat=(0.0,0.0,0.0, 0.0,0.0,0.0); 
  #purines ->#O6 sum_sqr, sum, n; #N7 sum_sqr, sum, n
  #pyrmimdines ->#O4 sum_sqr, sum, n; #O2 sum_sqr, sum, n
  
#=============================================================
#=============================================================
#  BEGIN MAIN PROGRAM
#-------------------------------------------------------------
   
while ($j < $r_time) { # operations
#while ($j < 10) { # testing
#   $frame_cnt = $frame_cnt + 1;
#=============================================================
# General File Operations
#-------------------------------------------------------------
#   $checkfile = $filein;
#   $filepdb = $base.$k;
#   $k=1;
   $filepdb=$base.'.'.$k; # open the file 

# see if the file exists 
   if (-e $filepdb) {
     $frame_cnt = $frame_cnt + 1;
     $attemptfile = 0;
     open(DATFILE, $filepdb) || die "Cannot open file $filepdb -- $!\n";
#  open(DATFILE, $filepdb) || attemptfile;
# case where the file skips a frame or two   
# attempt opening successive k files

   print "Processing output file ($filepdb)... \n";
   @lines = <DATFILE>;

   # basic statistics on occupancy
   $ions_inside = 0;
   $occu_cnt = 0; # reset ion count
   
   
   $atom_cnt = 0; # integer
   $set_cnt = 0; # reset the frame count for each file
   $atom_total=0; # resent total atoms

#=============================================================
# Phosphorous Atomic X,Y,Z & Holder variables for the ions
# ION X,Y,Z
#-------------------------------------------------------------
   $p1x=0.0;
   $p1y=0.0;
   $p1x=0.0;
 
   $p2x=0.0;
   $p2y=0.0;
   $p2x=0.0;

   $p3x=0.0;
   $p3y=0.0;
   $p3x=0.0;

   $p4x=0.0;
   $p4y=0.0;
   $p4x=0.0;
   
#check variables
   $p1c=0.0;
   $p2c=0.0;
   $p3c=0.0;
   $p4c=0.0;

   $m0x=0.0; # midpoint of cylinder
   $m0y=0.0; #
   $m0z=0.0; #

   $t0x=0.0; # end zero of the cylinder
   $t0y=0.0; 
   $t0z=0.0; 

   $t1x=0.0; # the other end of the cylinder
   $t1y=0.0; 
   $t1z=0.0;
#starts from zero indexing
#   open(MYFILE, '>'.$base.'out')|| die; # overwrite old files
#   close(MYFLIE);
##############################################################
# BEGIN LINE PROCESSING OF A FILE                            #
# Assumes that the main atoms (phosphorous etc) occur before #
# the atom identifications of the ions                       #
##############################################################
   @mid_cylinder = (0.0,0.0,0.0);
   @midp12 = (0.0,0.0,0.0);
   @midp34 = (0.0,0.0,0.0);
   $mid_made = 0.0; # assigned 1 after the final phosphorous is found
   foreach $line(@lines) {


   $n1x=0.0; # ion x
   $n1y=0.0; # ion y
   $n1z=0.0; # ion z
   @n1= (0.0,0.0,0.0);
   $d1c=0.0; # vertical to point X0 from the ray X21
   @r1c=(0.0,0.0,0.0,0.0); # radial distance from Xm to X0

# Search for the Phosphorus first
# Calculate the Cylinder Parameters
# Search for the Ions Next
# Calculate the line to the mid section of the cylider and
# distance from the mid section of the channel

       # throw out the first data after the "Processing Amber trajectory"

       @parsed_values = split(/\s+/, $line);
#       $atom_numid = trim($parsed_values[1])*1.0;
       # was set to $parsed_values[3]=~ /atoms/
        if ($parsed_values[0] =~ /TER/) { 
          #print "skip on TER \n"; 
        }
        elsif ($parsed_values[0] =~ m/ATOM/) {   # ensure that we are checking an atomic position
          
          $atom_numid = ($parsed_values[1])*1.0;

          if ($phos1 < 0 ) { # entering arbitrary cylinder mode when phos1 is negative
#              print "Arbitrary Cyilnder Mode \ \n"; # define a midpoint with the remaing three phosphorous id's
                                                    #x, y, z, the axis is along the z direction
#              print "Cylinder Center at x,y,z ->\ \ $phos2 \ \ , \ \ $phos3 \ \ , \ \  $phos4 \ \n";
              $p1x=$phos2 - $r_max;
              $p1y=$phos3 - $r_max;
              $p1z=$phos4 - $z_max;
              @p1=($p1x,$p1y,$p1z); # array
              $p1c=1.0;
              $p2x=$phos2 + $r_max;
              $p2y=$phos3 + $r_max;
              $p2z=$phos4 - $z_max;
              @p2=($p2x,$p2y,$p2z); # array
              $p2c=1.0;
              $p3x=$phos2 - $r_max;
              $p3y=$phos3 - $r_max;
              $p3z=$phos4 + $z_max;
              @p3=($p3x,$p3y,$p3z); # array
              $p3c=1.0;
              $p4x=$phos2 + $r_max;
              $p4y=$phos3 + $r_max;
              $p4z=$phos4 + $z_max;
              @p4=($p4x,$p4y,$p4z); # array
              $p4c=1.0;
              $mid_made = 1.0;
#              $atom_numid = ($parsed_values[1])*1.0;
          } else {
          #==============================================================
          # NORMAL OPERATIONS
          #--------------------------------------------------------------
#          $atom_numid = ($parsed_values[1])*1.0; 
#          print "read \ $parsed_values[1] \ \ $atom_numid \ \ $phos1 \n";
           if ($parsed_values[2]=~ m/P/ ) {
              if ($atom_numid eq $phos1) {
                 $p1x=($parsed_values[5])*1.0;
                 $p1y=($parsed_values[6])*1.0;
                 $p1z=($parsed_values[7])*1.0;
                 @p1=($p1x,$p1y,$p1z); # array 
                 $p1c = 1.0;
                 print "got 1 @p1 \n";
              } elsif ($atom_numid eq $phos2) {
                 $p2x=($parsed_values[5])*1.0;
                 $p2y=($parsed_values[6])*1.0;
                 $p2z=($parsed_values[7])*1.0;
                 @p2=($p2x,$p2y,$p2z); # array 
                 $p2c = 1.0;
                 print "got 2 @p2 \n";
              } elsif ($atom_numid eq $phos3){
                 $p3x=($parsed_values[5])*1.0;
                 $p3y=($parsed_values[6])*1.0;
                 $p3z=($parsed_values[7])*1.0;
                 @p3=($p3x,$p3y,$p3z); # array 
                 $p3c = 1.0;
                 print "got 3 @p3 \n";
              } elsif ($atom_numid eq $phos4){
                 $p4x=($parsed_values[5])*1.0;
                 $p4y=($parsed_values[6])*1.0;
                 $p4z=($parsed_values[7])*1.0; 
                 @p4=($p4x,$p4y,$p4z); # array 
                 $p4c = 1.0;
                 $mid_made = 1.0;
                 print "got 4 @p4 \n"; 
               }  # end phos search          
            #dump last count to the tally and append the datafile
          } # end of search for phosphorous
          
         } # end encapsulation of arbitrary cylinder mode

          # check if we have all the phosphorous coordinates
          # to begin calculating the cylinder dimensions

     
          #======================================
          # Create Midpoint and extrapolate
          # axis
          #--------------------------------------
          if ($p1c*$p2c*$p3c*$p4c*$mid_made gt 0.0) {
          # call up 
             @midp12=&get_midpt(@p1,@p2); # get middle of 1 2
             @midp34=&get_midpt(@p3,@p4); # get middle of 3 4
          # calculate the mid point of the channel
             @mid_cylinder=&get_midpt(@midp12,@midp34); # get middle of middle
             print "p12 \ @midp12 \ ,p34 \ @midp34 \ ,center-> \ @mid_cylinder \n";
             @axis_u_vec = get_u_vec(@midp12,@midp34); # get unit vector pointing parallel to the axis
             $mid_made = -1.0; # switch off once the mid point has been calculated
#             print "made midpoint \ @u_vec \n";
          }

          #=====================================
          # Populate full atom position, ID array
          # for distance and association matrices
          #-------------------------------------  
          if ($parsed_values[0] =~ m/TER/) { 
             print "skip on TER data cull \n"; } # null operation         
          elsif ($atom_numid < $ion_start) {
                 # array indices being at zero for perl
                 # record atom positions for distance measurements
                 $d_idx = $parsed_values[1]*1.0; #data index
                 $d_idx_nn1 = ($d_idx - 1)*5;

                 $atom_tab[$d_idx_nn1 + 0] =($parsed_values[5])*1.0; # x
                 $atom_tab[$d_idx_nn1 + 1] =($parsed_values[6])*1.0; # y
                 $atom_tab[$d_idx_nn1 + 2] =($parsed_values[7])*1.0; # z
                 # record residue type, unassigned is T which is in box zero
                 if    ($parsed_values[3]=~ m/G/) { $atom_tab[$d_idx_nn1 + 3] = 1; }
                 elsif ($parsed_values[3]=~ m/A/) { $atom_tab[$d_idx_nn1 + 3] = 2; }
                 elsif ($parsed_values[3]=~ m/C/) { $atom_tab[$d_idx_nn1 + 3] = 3; }
                 elsif ($parsed_values[3]=~ m/U/) { $atom_tab[$d_idx_nn1 + 3] = 4; }
                 
                 
                 $atom_tab[$d_idx_nn1 + 4] = $parsed_values[4]*1.0; # RESID Number

                 $res_max= $parsed_values[4]*1.0;   

                 # atom type array for statistics # was $d_idx_nn1
                 $atom_type_name[$d_idx]=trim($parsed_values[2]); # store the atom type

             #      if ($atom_type_name[$d_idx_nn1]=~ m/O1P/){
             #        print " got $atom_type_name[$d_idx_nn1] $d_idx_nn1 \n";
             #        exit;
             #      }
                  
               #   @atom_type_name=(0); # G:1,2-> O6,N7; A:3->N7
               #   @atom_type_stat=(0.0,0.0,0.0, 0.0,0.0,0.0);
               #   #purines ->#O6 sum_sqr, sum, n; #N7 sum_sqr, sum, n
               #   #pyrmimdines ->#O4 sum_sqr, sum, n; #O2 sum_sqr, sum, n
    
          }

          #=====================================
          # Perform Analysis
          #-------------------------------------
          # check if atom id is the ion id
          if ($mid_made lt 0.0) {
             if ($atom_numid > $ion_start) {
                if ($atom_numid < $ion_stop) {
          #===============================================
          # ions are sequential in their residue order
          #                   print "ion \ $atom_numid \n";
          # Determine ION Sequence ID
                   $ion_numid = $ion_numid + 1; 
          #===============================================
          # Channel Statistics/Quantitative localizations
          #-----------------------------------------------
                   $n1x=$parsed_values[5]*1.0; # assign positional coordinates
                   $n1y=$parsed_values[6]*1.0; # 
                   $n1z=$parsed_values[7]*1.0; #
                   @n1=($n1x,$n1y,$n1z); # assign 
                   @x_vec =cross_prod(@n1,@mid_cylinder,@midp34);
                   # get the projection of the ion ray onto the ray from midpoint cylinder axis to P34 anchor
                   $d1c = dot_prod(@n1,@mid_cylinder,@midp34); # projection along the ray midpt->34
 
                   @r1c = get_r_dist(@n1,@mid_cylinder); # ray to the midpoint, magnitude is in index 3
                   $dzc = get_z_dist($r1c[3],$x_vec[3]);
#                   print "@x_vec \ ->  \ $dzc  \n";
                   if($x_vec[3] < $r_max) {
                      if ($dzc < $z_max) {
                      

                     #  print "ion $atom_numid \ @x_vec[3] \ $r_max ->  \ $dzc \ $z_max -> \ $d1c \n";

                      $ion_total = $ion_total + 1; # increment total ion count
                      $occu_cnt = $occu_cnt + 1; # increment ion count for this frame

                      # tabulate the grid count
                      $zr_index = get_zr_index($dzc,@x_vec[3],$d1c); # get the zr grid index
                      $ion_grid_cnt[$zr_index] = $ion_grid_cnt[$zr_index] + 1; # increment grid index
                 #     $ion_frame = $ion_frame + 1; # number of ions within this frame

                     @zr_raw = get_zr_raw($zr_index);
                 #    print "$zr_index -> x, y \ $zr_raw[0] \ $zr_raw[1] \ $ion_grid_cnt[$zr_index] \n"; 
                 # 
                     open(MYFILE, '>>'.$base.'.'.'out')|| die;
                     print MYFILE "$frame_cnt\ $atom_numid \ $x_vec[3]\ $dzc\ $d1c\ $ion_total\n";
                     close (MYFILE);

                     #=================================
                     # Calculate the phi angle relative to the 
                     # Xprod of midCylinder->P1 and midCylinder->mid(P3,P4)
                     # counter clockwise direction
                     # looking into the p34 direction / end view
                     #                   
                     #                   |  p34 
                     #                   |  /
                     #                   | /
                     #            p2<----x---->p1
                     #                  /|_|
                     #                 / |
                     #(end view)     p21 |
                     #
                     # direction _|_ axis going to the right p1 is zer0 degrees/ 0 radians
                     # direction _|_ axis toward p2 is 180degrees/pi radians
                     #---------------------------------
                     # get the unit vector perp to the axis and anchor point 1
                     @phi_vec = cross_prod(@p1, @mid_cylinder, @midp34); # 
                     # get the new point on the axis for the ion perp 

                        if ( $d1c < 0){ # make z negative
                            $dzc = $dzc*(-1.0);
                         }

                     @new_z_point = move_pt(@mid_cylinder,$axis_u_vec[0],$axis_u_vec[1],$axis_u_vec[2],$dzc);

                     @new_perp_point = move_pt(@new_z_point,$phi_vec[4],$phi_vec[5],$phi_vec[6],1.0);

                     # magnitude should be equal to x_vec[3];
                     @r_perp_ion = get_r_dist(@n1, @new_z_point); # ray perp to the axis

                      #make distant point in space
                      @up_one_pt = move_pt(@new_z_point,$axis_u_vec[0],$axis_u_vec[1],$axis_u_vec[2],1.0);
                      # determin angle sign
                      @p1_uperp_vec = cross_prod(@up_one_pt, @new_z_point, @new_perp_point); #points in direction of p1 anchor

                      @p1_phi_point = move_pt(@new_z_point,$p1_uperp_vec[4],$p1_uperp_vec[5],$p1_uperp_vec[6],1.0); 

                      $d_phi_plus_minus_val = dot_prod(@n1, @new_z_point, @new_perp_point);
                      #   print "p1 phi point @p1_phi_point   \n"; # points toward p1 _|_ to axis
                      #   print "anchor phi point @new_z_point \n"; # on the cylinder axis
                      #   print "perp point @new_perp_point \n"; # points _|_ to axis toward negative phi
                      #   print "phi vec @phi_vec \n";

                      # determine angle from vector pointing along the cross of phos1 x midp34                   
                      $d_new_perp = dot_prod(@n1, @new_z_point, @p1_phi_point);
                      

                      $raw_ratio = $d_new_perp/$r_perp_ion[3];
                      $phi_value = acos($d_new_perp/$r_perp_ion[3]); # is an even function
                                            
     
                      if ($d_phi_plus_minus_val > 0) {
                         # print "is negative $d_phi_plus_minus_val  \n";
                         $phi_value = $phi_value*(-1.0);
                         # print "pi dot val $d_new_perp , radius val $r_perp_ion[3] \n";
                      }

                      #$m_verify_perp = dot_prod(@midp34,@new_z_point,@new_perp_point);
                    #  $m_verify_perp = dot_prod(@midp34,@new_z_point,@p1_phi_point);                      

                   #   print "phi is $phi_value for $raw_ratio check perp $m_verify_perp \n\n";
                      #print " $midp34[0], $midp34[1], $midp34[2] \n";
 
                    # if ($atom_numid > 726) {exit;}
 
                     # }                


                     #================================
                     # find highest occupancy to a chosen residue
                     #--------------------------------
                     #   $ion_num_idx = $atom_numid - $ion_start; # has already been adjusted by minus 1
                        
                        for ($e_idx = 0; $e_idx < $res_B; $e_idx++) { #clear the loop
                           $ion_is_here_flag[($ion_numid - 1)*$res_B + $e_idx + 1] = 0;
                        }
                     #=================================
                     # find the closest atom & residue type
                     #---------------------------------
                         $go_find_nearest = 1.0;
                         @this_atom = (0.0,0.0,0.0);
                         $last_type = 0;
                         $special_type = 0;
                         $special_dist = $ia_dist; #3.2; # max bonding distance, was 5
                         $special_residue_bin = 0;
                         $residue_box_bin = 0; 
                         $last_res_dist = $ia_dist; #3.2; #$r_max; #nominal value is 5 angstroms for any meaningful electrostatic interactions
               

                         #print "entering search for nearest -> last residue $res_max at $r_max \n";
                         for ($d_idx= 0; $d_idx < $ion_start; $d_idx++) {  #atom loop
                                $this_atom[0] = $atom_tab[0 + 5*$d_idx]; #x
                                $this_atom[1] = $atom_tab[1 + 5*$d_idx]; #y
                                $this_atom[2] = $atom_tab[2 + 5*$d_idx]; #z
                                @nearest_residue = get_r_dist(@n1,@this_atom);  # get the smallest value greater than zero
                                
                                #$dummy_check  = $atom_tab[3 + 5*$d_idx];
                                $dummy_check2 = 1*$atom_tab[4 + 5*$d_idx]; # actual residue id

                                # ===================
                                # occupancy near residue at least below 4 Angstroms, was hard-coded, now in the ion_is_here_dist variable
                                # -------------------                               
                                if ($nearest_residue[3] > 0 && $nearest_residue[3] < $ion_is_here_dist && $dummy_check2 > 0) {
                                   $k_idx = ($ion_numid - 1)*$res_B + $dummy_check2;

                                    
                                   $ion_is_here_flag[$k_idx]  = 1; # doesn't matter how often

                                   #if ($ion_numid < 2 && $dummy_check2 < 2 ) {
                                   #print "idx $k_idx got it for $ion_numid $res_B $dummy_check2 \n"; # okay for one loop in a set
                                   #exit; 
                                   #}

                                  # if ($ion_numid < 2) {
                                  # print "idx $k_idx for $ion_numid didn't get it $dummy_check2 \n";
                                  # exit;
                                  # }

                                }

                                # ===================
                                # cylinder occupancy data
                                # -------------------
                                #print "$d_idx, nearest is $nearest_residue[3], $dummy_check, $dummy_check2, @this_atom \n";
                                if ($nearest_residue[3] > 0 && $nearest_residue[3] < $last_res_dist && $dummy_check2 > 0) {
                                    $last_res_dist = 1.0*$nearest_residue[3] ; # take the lowest value
                                    $residue_type_bin = 1*$atom_tab[3 + 5*$d_idx]; # record the type TGACU-> 0,1,2,3,4 into the type box
                                    $residue_box_bin =  1*$atom_tab[4 + 5*$d_idx]; # record residue id number
                                    $go_find_nearest = 0.0; # found the nearest atom within the magic distance
                                     if ($residue_box_bin < 1) {
                                          print "read fault\n"; 
                                         exit;}


                                    #collect atom statistics on the special atoms O6, N7, O4, O2
                                    #print "$atom_type_name[$d_idx] , $residue_box_bin, $d_idx \n";
                                     if ($last_res_dist < $special_dist) {
                                       # print "$special_dist \n";
                                        $atom_nomen = $atom_type_name[$d_idx];
                                    if ($atom_type_name[$d_idx]=~ m/O6/)     {
                                       $special_type = $residue_type_bin;
                                       $last_type = 1; 
                                       $special_dist = $last_res_dist;
                                       $special_residue_bin=$residue_box_bin;
                                       print "O6 dist $special_dist \n";
                                    }elsif ($atom_type_name[$d_idx]=~ m/N7/)  {
                                       $special_type = $residue_type_bin;
                                       $last_type = 2;
                                       $special_dist = $last_res_dist;
                                       $special_residue_bin=$residue_box_bin;
                                       print " N7 dist $special_dist \n";
                                    } elsif ($atom_type_name[$d_idx]=~ m/'/)   {
                                        print "before $atom_type_name[$d_idx] at $last_res_dist \n"; 
                                        $residue_box_bin = -1 ; # don't count back bone interactions
                                         $ion_is_here_flag[$k_idx]  = -1; # take back the flag
                                    } else {
                                       if ($atom_type_name[$d_idx]=~ m/O4/)    {
                                        $special_type = $residue_type_bin;
                                        $last_type = 3;
                                        $special_dist = $last_res_dist;
                                        $special_residue_bin=$residue_box_bin;
                                        print "O4 dist $special_dist\n";
                                       }
                                       elsif ($atom_type_name[$d_idx]=~ m/O2/) {
                                        $special_type = $residue_type_bin;
                                        $last_type = 4;
                                        $special_dist = $last_res_dist;
                                        $special_residue_bin=$residue_box_bin;
                                        print "O2 dist $special_dist\n";
                                       }
                                       elsif ($atom_type_name[$d_idx]=~ m/P/) { # is a phosphate backbone screening
                                          print "phosphate $atom_type_name[$d_idx] \n";
                                          $residue_box_bin = -1 ; # don't count phosphate screening ions
                                           $ion_is_here_flag[$k_idx]  = -1; # take back the flag
                                       }
                                       else  {#$last_type = 0;
                                         print "atomic $d_idx $atom_type_name[$d_idx] at $last_res_dist \n"; 
                                       #  if ($atom_type_name[$d_idx]=~ m/ /) {
                                       #  exit; }
                                       } #don't count toward being special 
                                    } 
                                    }

                                }     
                             #  }
                         } # end of atom loop
                            
                            #=============================================
                            # Tabulate occupancy cumulative, highest, etc.
                            #---------------------------------------------
                            for ($e_idx = 0; $e_idx < $res_B; $e_idx++) { #residue loop

                               $f_idx = ($ion_numid - 1)*$res_B + $e_idx + 1; # mapping function
    
                               if ($ion_is_here_flag[$f_idx] > 0 ) {
                                  $ion_is_here_cnt[$f_idx] = $ion_is_here_cnt[$f_idx] + 1;

                               #   print " got it \n";
                               #   exit;

                               } else {

                                  if ($ion_is_here_cnt[$f_idx] > $ion_is_here_high[$f_idx]) {
                                     $ion_is_here_high[$f_idx] = $ion_is_here_cnt[$f_idx];
                                    #print "$ion_num_idx lost it $ion_is_here_high[$ion_num_idx] \n";
                                    #exit;
                                  }
    
                                  $ion_is_here_tot[$f_idx] = $ion_is_here_tot[$f_idx] + $ion_is_here_cnt[$f_idx]; # accumulate

                                  $ion_is_here_cnt[$f_idx] = 0; # clears the current count
                                }

                            } # end residue loop

                            #=====================
                            # atom type statistics
                            #----------------------

                            if ($go_find_nearest < 1.0 & $residue_box_bin > 0.0) {
                            $ions_inside = $ions_inside + 1;
                            $residue_box[$residue_box_bin] = $residue_box[$residue_box_bin] + 1; #increment frequency count
                            $residue_type[$residue_type_bin] = $residue_type[$residue_type_bin] + 1; # only five types frequency
                                                                                                     #T is zero box
         print "ion $ion_numid at $last_res_dist Ang $atom_nomen -> res $residue_box_bin type $residue_type_bin | $k ps \n";
         $atom_nomen = "";
                                if ($special_residue_bin > 0) {
                                   $idx_type = 0;
                                   $idx_type = ($special_residue_bin - 1)*7; # index by residue ID
                                   print "$last_type is atom type at $special_dist of $special_residue_bin of $special_type \n";
                                   if ($special_residue_bin != $residue_box_bin) {print "non special O6,N7,O4,O2 binding!\n"; }

                                if ($last_type < 2.0 && $last_type > 0.0) {   # O6 G only
                                   #$atom_type_stat=(0.0,0.0,0.0, 0.0,0.0,0.0); }
                                   $atom_type_stat[$idx_type + 0] = $atom_type_stat[$idx_type + 0] + $special_dist*$special_dist;
                                   $atom_type_stat[$idx_type + 1] = $atom_type_stat[$idx_type + 1] + $special_dist;
                                   $atom_type_stat[$idx_type + 2] = $atom_type_stat[$idx_type + 2] + 1;
                                   $atom_type_stat[$idx_type + 6] = $special_type;
                                               
                                } elsif ($last_type < 3.0) { # N7 A,G
                                   $atom_type_stat[$idx_type + 3] = $atom_type_stat[$idx_type + 3] + $special_dist*$special_dist;
                                   $atom_type_stat[$idx_type + 4] = $atom_type_stat[$idx_type + 4] + $special_dist;
                                   $atom_type_stat[$idx_type + 5] = $atom_type_stat[$idx_type + 5] + 1;
                                   $atom_type_stat[$idx_type + 6] = $special_type;
                                   print "put into N7 \n";
                                } elsif ($last_type < 4.0){ # O4 type U,T
                                   $atom_type_stat[$idx_type + 0] = $atom_type_stat[$idx_type + 0] + $special_dist*$special_dist;
                                   $atom_type_stat[$idx_type + 1] = $atom_type_stat[$idx_type + 1] + $special_dist;
                                   $atom_type_stat[$idx_type + 2] = $atom_type_stat[$idx_type + 2] + 1;
                                   $atom_type_stat[$idx_type + 6] = $special_type;

                                } elsif ($last_type < 5.0) { # O2 type C,U,T
                                   $atom_type_stat[$idx_type + 3] = $atom_type_stat[$idx_type + 3] + $special_dist*$special_dist;
                                   $atom_type_stat[$idx_type + 4] = $atom_type_stat[$idx_type + 4] + $special_dist;
                                   $atom_type_stat[$idx_type + 5] = $atom_type_stat[$idx_type + 5] + 1;
                                   $atom_type_stat[$idx_type + 6] = $special_type;

                                } else { 
                                  print "not interested in this atom type\n";
                                }
                               #purines ->#O6 sum_sqr, sum, n; #N7 sum_sqr, sum, n
                               #pyrmimdines ->#O4 sum_sqr, sum, n; #O2 sum_sqr, sum, n
                               }

                             }


                            #====================================                 
                            # z - phi tables
                            #------------------------------------
                            # $dzc is signed by now
                            $idx_z_phi = get_z_phi_idx($dzc, $phi_value); # increment phi by 0.0872664625997 radians or 5 degrees

                            # composite table
                            $z_phi_table[$idx_z_phi]=$z_phi_table[$idx_z_phi] + 1;

                            # stem 1 table z_phi_table of those near the helices
                            if ($residue_box_bin > 0 && $residue_box_bin < ($res_A + 1) ) {
                              $z_phi_tableA[$idx_z_phi]=$z_phi_tableA[$idx_z_phi] + 1;
                            } elsif ($residue_box_bin > $res_A && $residue_box_bin < ($res_B + 1) ) {
                            # stem 2 table z_phi_table
                              $z_phi_tableB[$idx_z_phi]=$z_phi_tableB[$idx_z_phi] + 1;
                            } else { print "not in the box\n";  }
                            #print "made z phi index, $idx_z_phi for $dzc, $phi_value\n";
                            #@dummy_z_phi = get_z_phi_raw($idx_z_phi);
                            
                            #print "raw return z->,$dummy_z_phi[0] ,phi-> $dummy_z_phi[1], cnt-> $z_phi_table[$idx_z_phi]  \n";
                            #exit;
                            #------------------------------------
                            # end of z-phi tabulations
                            #------------------------------------
                            #-------------------------------
                 
                       } # greater than z_max
                    } # greater than r_max
                 
                 #====================================
                 # bulk statistics
                 #------------------------------------
                 # remains in the annular region
                   if ($x_vec[3] > $br_min && $x_vec[3] < 40) { # was greater than
                     $bulk_cnt = $bulk_cnt + 1;
                     $bulk_total = $bulk_total + 1;
                     $ion_residency_cnt[$ion_numid] = $ion_residency_cnt[$ion_numid] + 1;

                   } elsif ($dzc > $bz_min && $dzc < 47) {
                     $bulk_cnt = $bulk_cnt + 1;
                     $bulk_total = $bulk_total + 1;
                     $ion_residency_cnt[$ion_numid] = $ion_residency_cnt[$ion_numid] + 1;
                           
                   } 
                  #==================================
                  # case that is in the RNA zone
                  # both radial distance and axial line are 
                  # less than the boxed min
                  # the bulk ion count is updated entering the cylinder
                  #----------------------------------
                   if ($x_vec[3] < $br_min || $x_vec[3] > 40) { # was less than
                      if ($dzc < $bz_min || $dzc > 47 ) {

                        if ($ion_residency_cnt[$ion_numid] > $ion_high_residency_cnt[$ion_numid]) {
                        $ion_high_residency_cnt[$ion_numid] = $ion_residency_cnt[$ion_numid];
                        }
                      $ion_residency_cnt[$ion_numid] = 0; # clear the residency count
                      }

                   }# elsif ($x_vec[3] > $br_min || $dzc > $bz_min) {#travels farther out from the annular region

                     #  if ($ion_residency_cnt[$ion_numid] > $ion_high_residency_cnt[$ion_numid]) {
                     #   $ion_high_residency_cnt[$ion_numid] = $ion_residency_cnt[$ion_numid];
                     #   }
                     # $ion_residency_cnt[$ion_numid] = 0; # clear the residency count
                     # }

                 # } 
                 #====================================              
                 #------------------------------------
                 # end bulk statistics
                 #------------------------------------

                 #====================================
                 # cylinder statistics
                 #------------------------------------
                 # when leaving the cylinder, the occupancy
                 # count is moved to the the high counter array

                 if ($x_vec[3] > $br_min ) { # was greater than

                        if ($ion_cylinder_residency_cnt[$ion_numid] > $ion_cylinder_high_residency_cnt[$ion_numid]) {
                        $ion_cylinder_high_residency_cnt[$ion_numid] = $ion_cylinder_residency_cnt[$ion_numid];
                        }

                        $ion_cylinder_residency_cnt[$ion_numid] = 0; # clear the residency count

                   } elsif ($dzc > $bz_min ) { # a larger cylinder of +/- 47 Angstroms

                        if ($ion_cylinder_residency_cnt[$ion_numid] > $ion_cylinder_high_residency_cnt[$ion_numid]) {
                        $ion_cylinder_high_residency_cnt[$ion_numid] = $ion_cylinder_residency_cnt[$ion_numid];
                        }

                        $ion_cylinder_residency_cnt[$ion_numid] = 0; # clear the residency count

                   }
                  #==================================
                  # case that it is in the RNA zone
                  # both radial distance and axial line are 
                  # less than the boxed min, occupancy is updated
                  #----------------------------------
                   if ($x_vec[3] < $br_min) { # both cases must be satisfied for the ion to be in the cylinder

                      if ($dzc < $bz_min) {
                      $cylinder_cnt = $cylinder_cnt + 1;
                      $cylinder_total = $cylinder_total + 1;
                      $ion_cylinder_residency_cnt[$ion_numid] = $ion_cylinder_residency_cnt[$ion_numid] + 1;
                      }

                   }
                 #====================================
                 #------------------------------------
                 # end cylinder statistics
                 #------------------------------------
                   if ($cylinder_cnt != $occu_cnt) {
                       print "$x_vec[3] \ \ $dzc \ \ $z_max \ \ $bz_min  \n";
                       print "count inequality \n";
                       print "bulk cnt \ \ $bulk_cnt \n";
                       print "occu cnt \ \ $occu_cnt \n";
                       print "cylinder_cnt \ \ $cylinder_cnt \n";
                       exit;
                   } 

                } # atom ion inner functions
              } # atom ion functions
            } # middle made
 
      } # end of search for atoms  in a file
   } # end of the foreach lines in a file
      $k=$k + $t_step; # increment time is 50ps, may change this to a command line variable
      $j=$j + 1;


      #clear ion sequence identifier
      $ion_high_numid = $ion_numid;
      $ion_numid = 0;

      #appending the output file with the results
#      open(MYFILE, '>>'.$base.'out')|| die;
#      print MYFILE "$frame_cnt\ $atom_cnt\ $atom_total\ $set_cnt\n";
#      close (MYFILE); 
       
       print "ions in channel  \ $occu_cnt \n";
      
       # add key statistical data
       $sqr_occu_sum = $sqr_occu_sum + $occu_cnt*$occu_cnt; # running sum of squares
       $occu_cnt = 0; # clear ion count

# bulk statistics
       print "ions in bulk \ $bulk_cnt \n";
       $sqr_bulk_sum = $sqr_bulk_sum + $bulk_cnt*$bulk_cnt; # running sum of squares
       $bulk_cnt = 0; # clear bulk ion count

# cylinder statistics
       print "ions in cylinder \ $cylinder_cnt \n";
       $sqr_cylinder_sum = $sqr_cylinder_sum + $cylinder_cnt*$cylinder_cnt;
       $cylinder_cnt = 0;   

# non screening ions
       print "ions inside \ $ions_inside \n";
       $ions_inside = 0;

      close(DATFILE);

  } else {

      if ($attemptfile < 99) {
#       $attemptfile = $attemptfile + 1;
          print "attempting \ $filepdb  \  $attemptfile \n";
#
# Attempts are made until the terminal loop is reached
#
          if ($attemptfile < 1) {
             $k=$k - $t_step; # increment by one to see if it is going up
             $k=$k + 1;
           } else {
             $k=$k + 1;
           }
         $attemptfile = $attemptfile + 1;

       } else { 
          report_stats;
          exit;
          }
   }
}  # end of file series processing and
# if it hasn't exited will exit now
report_stats;
