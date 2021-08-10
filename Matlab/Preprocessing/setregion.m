%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function reg = setregion(boxlo, boxhi, boxtilt)

if nargin<2
    error("Invalid input arguments");
end

if length(boxlo) ~= 3
    error("length of boxlo must be 3");
end    
if length(boxhi) ~= 3
    error("length of boxhi must be 3");
end    

reg.boxlo = boxlo;
reg.boxhi = boxhi;

if nargin>=3 && ~isempty(boxtilt)
    if length(boxtilt) ~= 3
        error("length of boxtilt must be 3");
    end        
    reg.boxtilt = boxtilt;
    reg.triclinic = 1;
else
    reg.boxtilt = [0 0 0];
    reg.triclinic = 0;
end

end


% struct regionstruct {       
%     int style;
%     int interior;                     // 1 for interior, 0 for exterior
%     int scaleflag;                    // 1 for lattice, 0 for box  
%     int bboxflag;                // 1 if bounding box is computable
%     int varshape;                // 1 if region shape changes over time
%     int dynamic;                 // 1 if position/orient changes over time
%     int moveflag, rotateflag;    // 1 if position/orientation changes
%     int openflag;                // 1 if any face is open
%     
%     dstype *r;                   // distance between particle & surf, r > 0.0
%     dstype *delx, *dely, *delz;    // vector from surface pt to particle
%     dstype *radius;              // curvature of region at contact point
%     int *iwall;                  // unique id of wall for storing shear history    
%     int *openfaces;//[6];           // flags for which faces are open
%     int *tri;//[12*3];    // 3 corner pts of 12 triangles (2 per face)
%     
%     dstype *face;//[6*3];   // unit normals of 6 prism faces
%     dstype *corners;//[8*3];  // 8 corner pts of prism     
%     dstype *faces;//[6*4*3];   // 4 corner pts of 6 prism faces   
%     dstype *lo;//[3];
%     dstype *hi;//[3];
%     dstype *tilt;//[3];
%     dstype *clo;//[3]; // opposite corners of prism
%     dstype *chi;//[3]; // opposite corners of prism   
%     dstype *scale;//[3];
%     dstype *extent_lo;//[3];
%     dstype *extent_hi;//[3];
%     dstype *a;//[3];  // edge vectors of region
%     dstype *b;//[3];  // edge vectors of region
%     dstype *c;//[3];  // edge vectors of region
%     dstype *h;//[6];
%     dstype *h_inv;//[6];    
%     dstype *dx;//[3];    // current displacement 
%     dstype *v;//[3];                 // translational velocity
%     dstype *rpoint;//[3];            // current origin of rotation axis
%     dstype *omega;//[3];             // angular velocity    
%     dstype *xcenter;//[3];    // translated/rotated center of cylinder/sphere (only used if varshape)
%     dstype *prev;//[5];       // stores displacement (X3), angle and if
%     dstype theta; // orientation
%     dstype rprev; // speed of time-dependent radius, if applicable
%     

