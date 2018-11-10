#!/usr/local/bin/perl

#
# By Latsavongsakda Sethaphong
#    Yingling Lab
#    NCSU 

# The formulation and general algorithms are taken from 
# Alexandra Patriksson and David van der Spoel
# A temperature pedictor for parallel tempering simulations
# Physical Chemistry Chemical Physics 2008
# pp. 2073-2077.

# For efficient exchange between replicants, each replica
# should spend equal amounts of time at each temperature
# in order to facilitate structure interchanges from the
# highest temperatures to the lowest temperatures.
# The requirement is such that the acceptance probability
# between all adjacent pairs need to be equal.
# 
# the C parametrization is somewhat different from
# what is given in Patriksson & van der Spoel, but is
# within bounds for predicting likely temperatures
# for a given/desired exchange probability.
#
# The only important terms are $A1, the specific heat
# for water, $B1, the specific heat of the constituent  
# of interest.
#
# Version 1.0 30 July, 2009
#
##########################################################
use Math::Trig;

sub calc_mu_cv;
sub calc_Ndf;
sub calc_pmu_1_mu_2;
sub calc_mu_12;
sub calc_sigma_12;
sub calc_erf; # calculate the guassian error function
sub factorial; # utility subroutine
sub calc_dir; # step and directional routine
sub tabulate_temp; # main routine

sub calc_mucv {
 my ($Tnow) = @_;
 my $mu_cv;
 $mu_cv = 0.5*$kB*$Tnow * ($myNc + $myNv);

 return $mu_cv;

}

sub calc_Ndf {
my ($myNw, $myNp, $myNc, $myNv) = @_; # Take in the inputs
my $myNdf = 0.0;

$myNdf = 9*$myNw + 3*$myNp - $myNc - $myNv;
return $myNdf;
}

sub calc_pmu1_mu2 {
#transition probability
   my ($mu12, $sigma12, $myC) = @_;

   $mu12_over_sigma12_root2 = (((-1.0)*$mu12)/$sigma12)/sqrt(2);
   
   $Cmu12_plus_Csig12o2 = $myC*$mu12 + ($myC*$myC*$sigma12*$sigma12)*0.5;

   $mu12_plus_Csigma12_sqrd_over_sigma12_root2 = ($mu12 + $myC*$sigma12*$sigma12)/($sigma12*(sqrt(2)));

#   print "calling erf once \n";
#   $first_term  = 0.5*(1.0 + calc_erf($mu12_over_sigma12_root2, 100));
   $first_term = 0.5*(1.0 + calc_erf($mu12_over_sigma12_root2, 150) );

   $middle_term = 0.5*(exp($Cmu12_plus_Csig12o2)); # middle term

#   print "calling erf twice \n";
   $last_term   = (1.0 + calc_erf($mu12_plus_Csigma12_sqrd_over_sigma12_root2, 150) );

   $prob_trans_12 = $first_term + $middle_term*$last_term;
#   print "terms $first_term, $middle_term, $last_term \n";
#   print "inputs $mu12_over_sigma12_root2, $Cmu12_plus_Csig12o2, $mu12_plus_Csigma12_sqrd_over_sigma12_root2\n";
return $prob_trans_12;
}

sub calc_mu_12{

  my ($myT1, $myT2, $myNw, $myNp, $myNc, $myNv) = @_;
  
  $mu_12 = ($myT2 - $myT1)*($A1*$myNw + $B1*$myNp - 0.5*$kB*($myNc + $myNv)); # take the independent value of C
 
#  print "value of mu 12 $mu_12, $myT2, $myT1, $myNw, $myNp \n";
  return $mu_12;
}

sub calc_sigma_12{
  my ($myNdf, $myT1, $myT2) = @_;
  $temp_sum = 0.0;
  $temp_sum = $myNdf*( ($D1**2)*($myT1**2 + $myT2**2) + 2*$D1*$D0*($myT2 + $myT1) + 2*($D0**2) );
  
  $sigma_12 = 0.0;
  $sigma_12 = sqrt($temp_sum);
#  print "value of sigma 12 $sigma_12\n";
  return $sigma_12;
}

sub calc_erf_cap_pi {
 
  my ($myz_square,$myn_sub) = @_; # 
 # my $n_sub_l1 = 0;
  my $aval = 0;

  #print "$myn_sub is submitted to cap pi \n";

  if ($myn_sub < 1) {
    $aval = 1; # for the zero case
#    print " $aval is last cap pi\n";
    return $aval;
  } else  {
#    $n_sub_l1 = $n_sub - 1;
    $N_fac = factorial($myn_sub);
    $aval = ((-1.0)**$myn_sub)*($myz_square**$myn_sub)/$N_fac;
#    print " in cap_pi $aval for $myn_sub, $N_fac \n";
    #$aval= $myz_square*(-1.0)*(1/$n_sub)*calc_erf_cap_pi($myz_square, $n_sub_l1); # recursive call
    return $aval;
  }

}

sub calc_erf {
#
# Returns the gaussian error function
# The external call should scale by 2/sqrt(pi)
#
   my ($z_raw, $max_n) = @_; 
   my $cap_pi_val = 0;
   $pi = 3.1415926535897932384626433832795;
   $z_square = $z_raw*$z_raw;
   $z_over_odd_n = 0;
   my $n_less_1 = $max_n - 1;

#   if ($z_raw < 1) {
#    print "error on input in erf \n";
#    return 0;
#   }

   if ($max_n < 1) {
     #$cap_pi_val = calc_erf_cap_pi($z_square,0);
     $z_over_odd_n = $z_raw*2/sqrt($pi);
     #$z_over_odd_n = ($z_raw/(2*$max_n + 1))*$cap_pi_val; # for n= 0, k=1
#     print "in erf last $z_over_odd_n \n";
     return $z_over_odd_n;
   } else {
     $cap_pi_val = calc_erf_cap_pi($z_square, $max_n);
     $z_over_odd_n = ($z_raw/(2*$max_n + 1))*$cap_pi_val*2/(sqrt($pi)) + calc_erf($z_raw, $n_less_1); # recursive call
#     print "in erf $z_over_odd_n, $z_raw, $cap_pi_val, $max_n \n";
     return $z_over_odd_n;
   }
}


sub calc_C {
#
# Patriksson and van der Spoel give C as
# 1/kBT1 - 1/kBT2, which comes out to be a very large number
# consequenty, it doesn't work. An email requesting clarification
# has been sent. No reply yet
#
  my($myT1, $myT2) = @_;
  my $myC_val = 0.0;

  # power of -0.30 works well
#   $myC_val =  -( ( 1/($kB*$myT1) -  1/($kB*$myT2))**(-0.298) ); # this needs to be negative

   $myC_val = ( 1/(8.314472*$myT2) - 1/(8.314472*$myT1)); # gas constant
#  print "C is  $myC_val  \n";
  return $myC_val;
}

sub factorial {
  my($num) = @_;
  if($num ==  1) { 
    return 1; # stops at 1, factorial does not multiply zero
  } else {
    return $num*factorial($num - 1); # recursive call
  }
}

sub calc_dir {
#
# If the probability is too high, increment Temperature separation
# If the probability is too low, decrement Temperature separation

   my ($myPdes, $myPpred, $step) = @_;
   $myPdes = $pdes;
   my $new_step = (($myPpred - $myPdes)/($myPdes))*$step; 

return $new_step;
}


# iteration limit for convergence
sub tabulate_temp {

 my ($start_Temp,$stop_Temp) = @_;
 my $iter_limit = 200;
 my $i = 1;
 my $first_temp = 0;
 my $k_step = 0.25;
 my $new_step = 1.0;
 my $mu_12 = 1.0;
 my $sigma_12 = 1.0;
 my $C_val = 1.0;
 my $P12 = 0.0;
 my $tolerance =0.0; 
  
   $first_Temp = $start_Temp;
   $next_Temp = $start_Temp;
   $new_step = $k_step;

   $Ndf = calc_Ndf($num_wat, $num_p, $num_c, $num_v); #water, atoms of nucleic acid, constraints, virtual sites

   open(MYFILE, '>'.$f_out.'.'.'dat');
   close(MYFILE);

   open(MYFILE, '>>'.$f_out.'.'.'dat')||die;
   printf MYFILE "%-5.1f", $start_Temp;
   print MYFILE "\n";

   print "$Ndf is the number of dregrees of freedom.\n";

 while ($next_Temp < $stop_Temp) {
   $i=1;
   while ($i < $iter_limit) {
      $next_Temp = $next_Temp + $new_step;

      $C_val = calc_C($first_Temp, $next_Temp); # the starting temperature values
      $mu_12 = calc_mu_12($first_Temp, $next_Temp, $num_wat, $num_p,$num_v, $num_c);
      $sigma_12 = calc_sigma_12($Ndf, $first_Temp, $next_Temp);
      $P12 = calc_pmu1_mu2($mu_12, $sigma_12,$C_val);

      if ($P12 < $pdes) { # overshoot first then, draw back
         $new_step = calc_dir($pdes,$P12,$k_step);
#         print "over shot $P12 < $pdes \n";
      }# else { print "under shot $P12 > $pdes \n";}
      $i = $i + 1;

      $tolerance = $new_step/$k_step;

      if ($i == $iter_limit) {
      print " reached iteration limit\n";
      print "tol $tolerance for P($first_Temp,$next_Temp) = $P12 \n";
      exit;
      }
      
      if (abs($tolerance) < abs($num_tol)) {
       print "met tol $tolerance for P($first_Temp,$next_Temp) = $P12 \n";
       printf MYFILE "%-5.1f", $next_Temp;
       print MYFILE "\n";
       $i = $iter_limit;
      # exit;
      }

    
    } #inner while loop
#    return $fin_temp;
 
   $first_Temp = $next_Temp;
   $new_step = $k_step;

  } # outer while loop
  close(MYFILE);
  print " completed \n";
  return;
}


####################################
#
# BEGIN MAIN ROUTINES 
#
#
###################################

if ($#ARGV < 0) {
   print " Incorrect usage...\n";
   exit;
}

  $j=1;
  $f_out   = $ARGV[0]; # file name output (e.g. temp.dat)
  $pdes    = $ARGV[1]*1.0; # desired probability [0.0,1.0]
  $templ   = $ARGV[2]*1.0; # lower temperature (kelvin)
  $tempu   = $ARGV[3]*1.0; # upper temperature (kelvin)
  $num_wat = $ARGV[4]*1.0; # number of waters
  $num_p   = $ARGV[5]*1.0; # number of particles
  $num_c   = $ARGV[6]*1.0; # number of constraints
  $num_v   = $ARGV[7]*1.0; # number of virtual sites
  $num_tol = $ARGV[8]*1.0; # relative convergence tolerance (small like 0.001)
  $mode    = $ARGV[9]*1.0; # Protein(3), DNA(1), RNA(2), default is Protein
# basic equations of relevance
# mu_12(T)= (T2-T1){ A1.Nw +B1.Np - 0.5(kb)(Nc + Nv)}
# sigma_12(T)=sqrt{ Ndf{D^2 [(T2)^2 + (T1)^2] + 2D1.D0.[T2 + T1] + 2(D0)^2 }
# Ndf = 9.Nw + 3.Np - (Nc + Nv)
# mu(T) = (A0 + A1.T)Nw + (B0 + B1.T)Np - 0.5(kb)T(Nc + Nv)
# sigma(T) = (D0 + D1.T)sqr(Ndf)
# pu12 = {1/[sigma_12*sqr(2*pi)]} * exp[{mu(T)-mu_12(T)}/(2*sigma_12^2)
# 
# Fitting for the temperature intervals are approximate
#
#
#  Everything is in Joules not k-Joules
#
#   $A0 = -59.22*1000; # kJ/mol
#   $A1 =  75.94; # J/mol-K
#   $B0 = -22.84*1000; # kJ/mol
#   $B1 =  13.47; # J/mol-K
#   $D0 =  1.169*1000; # kJ/mol
#   $D1 =  2.976; # J/mol-K

# Boltzmann's constant
$kB = 1.3806503*(10**(-23)); # m^2 kg s^-2 K^-2 (J/K-1)

########################
# Set Linear Constants #
########################

 if ($mode > 0) {
   if ($mode = 1) { # DNA
   $A0 = -59.22*1000; # kJ/mol
   $A1 =  75.94; # J/mol-K
   $B0 = -22.84*1000; # kJ/mol
   $B1 =  13.47; # J/mol-K
   $D0 =  1.169*1000; # kJ/mol
   $D1 =  2.976; # J/mol-K

   } elsif ($mode = 2) { # RNA


#############################
   $A0 = -59.22*1000; # kJ/mol
   $A1 =  75.94; # J/mol-K
   $B0 = -22.84*1000; # kJ/mol
   $B1 =  13.47; # J/mol-K
   $D0 =  1.169*1000; # kJ/mol
   $D1 =  2.976; # J/mol-K

# Protein by default        #
############################
   } else { # from A. Patriksson & D. van der Spoel

   $A0 = -59.22*1000; # kJ/mol
   $A1 =  75.94; # J/mol-K
   $B0 = -22.84*1000; # kJ/mol
   $B1 =  13.47; # J/mol-K
   $D0 =  1.169*1000; # kJ/mol
   $D1 =  2.976; # J/mol-K

   }
 } else { # from A. Patriksson & D. van der Spoel

   $A0 = -59.22*1000; # kJ/mol
   $A1 =  75.94; # J/mol-K
   $B0 = -22.84*1000; # kJ/mol
   $B1 =  13.47; # J/mol-K
   $D0 =  1.169*1000; # kJ/mol
   $D1 =  2.976; # J/mol-K


 }

 print " begin num tol $num_tol \n";
 tabulate_temp($templ,$tempu);
# $test_erf = calc_erf(3.5, 100);
#  print "estimated $test_erf \n";
