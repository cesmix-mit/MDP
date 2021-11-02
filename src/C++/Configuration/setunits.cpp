/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

void setunits(commonstruct &common, int style)
{
  if (style == 0) { // lj
    common.boltz = 1.0;
    common.hplanck = 1.0;
    common.mvv2e = 1.0;
    common.ftm2v = 1.0;
    common.mv2d = 1.0;
    common.nktv2p = 1.0;
    common.qqr2e = 1.0;
    common.qe2f = 1.0;
    common.vxmu2f = 1.0;
    common.xxt2kmu = 1.0;
    common.e_mass = 0.0;    // not yet set
    common.hhmrr2e = 0.0;
    common.mvh2r = 0.0;
    common.angstrom = 1.0;
    common.femtosecond = 1.0;
    common.qelectron = 1.0;
    common.dt = 0.005;
    common.skin = 0.3;

  } else if (style == 1) { // real
    common.boltz = 0.0019872067;
    common.hplanck = 95.306976368;
    common.mvv2e = 48.88821291 * 48.88821291;
    common.ftm2v = 1.0 / 48.88821291 / 48.88821291;
    common.mv2d = 1.0 / 0.602214129;
    common.nktv2p = 68568.415;
    common.qqr2e = 332.06371;     // see also common.qqr2d_lammps_real
    common.qe2f = 23.060549;
    common.vxmu2f = 1.4393264316e4;
    common.xxt2kmu = 0.1;
    common.e_mass = 1.0/1836.1527556560675;
    common.hhmrr2e = 0.0957018663603261;
    common.mvh2r = 1.5339009481951;
    common.angstrom = 1.0;
    common.femtosecond = 1.0;
    common.qelectron = 1.0;

    common.dt = 1.0;
    common.skin = 2.0;

  } else if (style == 2) { // metal
    common.boltz = 8.617343e-5;
    common.hplanck = 4.135667403e-3;
    common.mvv2e = 1.0364269e-4;
    common.ftm2v = 1.0 / 1.0364269e-4;
    common.mv2d = 1.0 / 0.602214129;
    common.nktv2p = 1.6021765e6;
    common.qqr2e = 14.399645;
    common.qe2f = 1.0;
    common.vxmu2f = 0.6241509647;
    common.xxt2kmu = 1.0e-4;
    common.e_mass = 0.0;    // not yet set
    common.hhmrr2e = 0.0;
    common.mvh2r = 0.0;
    common.angstrom = 1.0;
    common.femtosecond = 1.0e-3;
    common.qelectron = 1.0;

    common.dt = 0.001;
    common.skin = 2.0;

  } else if (style == 3) { // si
    common.boltz = 1.3806504e-23;
    common.hplanck = 6.62606896e-34;
    common.mvv2e = 1.0;
    common.ftm2v = 1.0;
    common.mv2d = 1.0;
    common.nktv2p = 1.0;
    common.qqr2e = 8.9876e9;
    common.qe2f = 1.0;
    common.vxmu2f = 1.0;
    common.xxt2kmu = 1.0;
    common.e_mass = 0.0;    // not yet set
    common.hhmrr2e = 0.0;
    common.mvh2r = 0.0;
    common.angstrom = 1.0e-10;
    common.femtosecond = 1.0e-15;
    common.qelectron = 1.6021765e-19;
    common.dt = 1.0e-8;
    common.skin = 0.001;
  } else if (style == 4) { // cgs
    common.boltz = 1.3806504e-16;
    common.hplanck = 6.62606896e-27;
    common.mvv2e = 1.0;
    common.ftm2v = 1.0;
    common.mv2d = 1.0;
    common.nktv2p = 1.0;
    common.qqr2e = 1.0;
    common.qe2f = 1.0;
    common.vxmu2f = 1.0;
    common.xxt2kmu = 1.0;
    common.e_mass = 0.0;    // not yet set
    common.hhmrr2e = 0.0;
    common.mvh2r = 0.0;
    common.angstrom = 1.0e-8;
    common.femtosecond = 1.0e-15;
    common.qelectron = 4.8032044e-10;
    common.dt = 1.0e-8;
    common.skin = 0.1;

  } else if (style == 5) { // electron
    common.boltz = 3.16681534e-6;
    common.hplanck = 0.1519829846;
    common.mvv2e = 1.06657236;
    common.ftm2v = 0.937582899;
    common.mv2d = 1.0;
    common.nktv2p = 2.94210108e13;
    common.qqr2e = 1.0;
    common.qe2f = 1.94469051e-10;
    common.vxmu2f = 3.39893149e1;
    common.xxt2kmu = 3.13796367e-2;
    common.e_mass = 0.0;    // not yet set
    common.hhmrr2e = 0.0;
    common.mvh2r = 0.0;
    common.angstrom = 1.88972612;
    common.femtosecond = 1.0;
    common.qelectron = 1.0;

    common.dt = 0.001;
    common.skin = 2.0;

  } else if (style == 6) { // micro
    common.boltz = 1.3806504e-8;
    common.hplanck = 6.62606896e-13;
    common.mvv2e = 1.0;
    common.ftm2v = 1.0;
    common.mv2d = 1.0;
    common.nktv2p = 1.0;
    common.qqr2e = 8.987556e6;
    common.qe2f = 1.0;
    common.vxmu2f = 1.0;
    common.xxt2kmu = 1.0;
    common.e_mass = 0.0;    // not yet set
    common.hhmrr2e = 0.0;
    common.mvh2r = 0.0;
    common.angstrom = 1.0e-4;
    common.femtosecond = 1.0e-9;
    common.qelectron = 1.6021765e-7;
    common.dt = 2.0;
    common.skin = 0.1;

  } else if (style == 7) { // nano
    common.boltz = 0.013806504;
    common.hplanck = 6.62606896e-4;
    common.mvv2e = 1.0;
    common.ftm2v = 1.0;
    common.mv2d = 1.0;
    common.nktv2p = 1.0;
    common.qqr2e = 230.7078669;
    common.qe2f = 1.0;
    common.vxmu2f = 1.0;
    common.xxt2kmu = 1.0;
    common.e_mass = 0.0;    // not yet set
    common.hhmrr2e = 0.0;
    common.mvh2r = 0.0;
    common.angstrom = 1.0e-1;
    common.femtosecond = 1.0e-6;
    common.qelectron = 1.0;
    common.dt = 0.00045;
    common.skin = 0.1;

  } else error("Illegal units command");
}

