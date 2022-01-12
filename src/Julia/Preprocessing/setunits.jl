mutable struct UnitsStruct

    style::String; 
    boltz::Float64;
    hplanck::Float64;
    mvv2e::Float64;
    ftm2v::Float64;
    mv2d::Float64;
    nktv2p::Float64;
    qqr2e::Float64;
    qe2f::Float64;
    vxmu2f::Float64;
    xxt2kmu::Float64;
    e_mass::Float64;    # not yet set
    hhmrr2e::Float64;
    mvh2r::Float64;
    angstrom::Float64;
    femtosecond::Float64;
    qelectron::Float64;
    dt::Float64;
    skin::Float64;

    UnitsStruct() = new();
end

function setunits(style)

style = lowercase(style)
units = UnitsStruct();    
units.style = style 

if (style == "lj")  # lj
    units.boltz = 1.0;
    units.hplanck = 1.0;
    units.mvv2e = 1.0;
    units.ftm2v = 1.0;
    units.mv2d = 1.0;
    units.nktv2p = 1.0;
    units.qqr2e = 1.0;
    units.qe2f = 1.0;
    units.vxmu2f = 1.0;
    units.xxt2kmu = 1.0;
    units.e_mass = 0.0;    # not yet set
    units.hhmrr2e = 0.0;
    units.mvh2r = 0.0;
    units.angstrom = 1.0;
    units.femtosecond = 1.0;
    units.qelectron = 1.0;
    units.dt = 0.005;
    units.skin = 0.3;

elseif (style == "real")  # real
    units.boltz = 0.0019872067;
    units.hplanck = 95.306976368;
    units.mvv2e = 48.88821291 * 48.88821291;
    units.ftm2v = 1.0 / 48.88821291 / 48.88821291;
    units.mv2d = 1.0 / 0.602214129;
    units.nktv2p = 68568.415;
    units.qqr2e = 332.06371;     # see also units.qqr2d_lammps_real
    units.qe2f = 23.060549;
    units.vxmu2f = 1.4393264316e4;
    units.xxt2kmu = 0.1;
    units.e_mass = 1.0/1836.1527556560675;
    units.hhmrr2e = 0.0957018663603261;
    units.mvh2r = 1.5339009481951;
    units.angstrom = 1.0;
    units.femtosecond = 1.0;
    units.qelectron = 1.0;

    units.dt = 1.0;
    units.skin = 2.0;

elseif (style == "metal")  # metal
    units.boltz = 8.617343e-5;
    units.hplanck = 4.135667403e-3;
    units.mvv2e = 1.0364269e-4;
    units.ftm2v = 1.0 / 1.0364269e-4;
    units.mv2d = 1.0 / 0.602214129;
    units.nktv2p = 1.6021765e6;
    units.qqr2e = 14.399645;
    units.qe2f = 1.0;
    units.vxmu2f = 0.6241509647;
    units.xxt2kmu = 1.0e-4;
    units.e_mass = 0.0;    # not yet set
    units.hhmrr2e = 0.0;
    units.mvh2r = 0.0;
    units.angstrom = 1.0;
    units.femtosecond = 1.0e-3;
    units.qelectron = 1.0;

    units.dt = 0.001;
    units.skin = 2.0;

elseif (style == "si")  # si
    units.boltz = 1.3806504e-23;
    units.hplanck = 6.62606896e-34;
    units.mvv2e = 1.0;
    units.ftm2v = 1.0;
    units.mv2d = 1.0;
    units.nktv2p = 1.0;
    units.qqr2e = 8.9876e9;
    units.qe2f = 1.0;
    units.vxmu2f = 1.0;
    units.xxt2kmu = 1.0;
    units.e_mass = 0.0;    # not yet set
    units.hhmrr2e = 0.0;
    units.mvh2r = 0.0;
    units.angstrom = 1.0e-10;
    units.femtosecond = 1.0e-15;
    units.qelectron = 1.6021765e-19;
    units.dt = 1.0e-8;
    units.skin = 0.001;
elseif (style == "cgs")  # cgs
    units.boltz = 1.3806504e-16;
    units.hplanck = 6.62606896e-27;
    units.mvv2e = 1.0;
    units.ftm2v = 1.0;
    units.mv2d = 1.0;
    units.nktv2p = 1.0;
    units.qqr2e = 1.0;
    units.qe2f = 1.0;
    units.vxmu2f = 1.0;
    units.xxt2kmu = 1.0;
    units.e_mass = 0.0;    # not yet set
    units.hhmrr2e = 0.0;
    units.mvh2r = 0.0;
    units.angstrom = 1.0e-8;
    units.femtosecond = 1.0e-15;
    units.qelectron = 4.8032044e-10;
    units.dt = 1.0e-8;
    units.skin = 0.1;

elseif (style == "electron")  # electron
    units.boltz = 3.16681534e-6;
    units.hplanck = 0.1519829846;
    units.mvv2e = 1.06657236;
    units.ftm2v = 0.937582899;
    units.mv2d = 1.0;
    units.nktv2p = 2.94210108e13;
    units.qqr2e = 1.0;
    units.qe2f = 1.94469051e-10;
    units.vxmu2f = 3.39893149e1;
    units.xxt2kmu = 3.13796367e-2;
    units.e_mass = 0.0;    # not yet set
    units.hhmrr2e = 0.0;
    units.mvh2r = 0.0;
    units.angstrom = 1.88972612;
    units.femtosecond = 1.0;
    units.qelectron = 1.0;

    units.dt = 0.001;
    units.skin = 2.0;

elseif (style == "micro")  # micro
    units.boltz = 1.3806504e-8;
    units.hplanck = 6.62606896e-13;
    units.mvv2e = 1.0;
    units.ftm2v = 1.0;
    units.mv2d = 1.0;
    units.nktv2p = 1.0;
    units.qqr2e = 8.987556e6;
    units.qe2f = 1.0;
    units.vxmu2f = 1.0;
    units.xxt2kmu = 1.0;
    units.e_mass = 0.0;    # not yet set
    units.hhmrr2e = 0.0;
    units.mvh2r = 0.0;
    units.angstrom = 1.0e-4;
    units.femtosecond = 1.0e-9;
    units.qelectron = 1.6021765e-7;
    units.dt = 2.0;
    units.skin = 0.1;

elseif (style == "nano")  # nano
    units.boltz = 0.013806504;
    units.hplanck = 6.62606896e-4;
    units.mvv2e = 1.0;
    units.ftm2v = 1.0;
    units.mv2d = 1.0;
    units.nktv2p = 1.0;
    units.qqr2e = 230.7078669;
    units.qe2f = 1.0;
    units.vxmu2f = 1.0;
    units.xxt2kmu = 1.0;
    units.e_mass = 0.0;    # not yet set
    units.hhmrr2e = 0.0;
    units.mvh2r = 0.0;
    units.angstrom = 1.0e-1;
    units.femtosecond = 1.0e-6;
    units.qelectron = 1.0;
    units.dt = 0.00045;
    units.skin = 0.1;
else 
    error("Illegal units command");

end

return units 

end

function getunits(units)
    boltz = units.boltz
    hplanck = units.hplanck
    mvv2e = units.mvv2e
    ftm2v = units.ftm2v
    mv2d = units.mv2d
    nktv2p = units.nktv2p
    qqr2e = units.qqr2e
    qe2f = units.qe2f
    vxmu2f = units.vxmu2f
    xxt2kmu = units.xxt2kmu
    e_mass = units.e_mass    # not yet set
    hhmrr2e = units.hhmrr2e
    mvh2r = units.mvh2r
    angstrom = units.angstrom
    femtosecond = units.femtosecond
    qelectron = units.qelectron
    dt = units.dt
    skin = units.skin
end