//+
SetFactory("OpenCASCADE");
//+geometry
l=100;
n=500;
a_max=-200;
a_min=-700;
el=1000; //
f=25; //+dx of fracture
//+
Box(1) = {a_min, -0.5*el, -0.5*el, n, el, el};
//+
Box(2) = {a_min+2*l, -0.5*el, -0.5*el, l, el, el};
//+
Box(3) = {a_min, -0.5*el, -52.5, n, el, 105};
//+
Box(4) = {a_min, -52.5, -0.5*el, n, 105, el};
//+
BooleanFragments{ Volume{1,2,3,4}; Delete; }{}
//+

Transfinite Curve{137,115,79,76} = 22;
//+
Transfinite Curve{126,92,48,45} = 22;
//+
Transfinite Curve{105,64,25,21} = 22;
//+
Transfinite Curve{108,67,28,23} = 22;
//+
Transfinite Curve{87,90,124,141} = 22;
//+
Transfinite Curve{49,46,93,127} = 22;
//+
Transfinite Curve{15,20,54,98} = 22;
//+
Transfinite Curve{13,16,50,94} = 22;

//+
Transfinite Curve{133,111,83,82,114} = 21;
Transfinite Curve{122,99,91,63,52,36,51,31} = 21;
Transfinite Curve{101,65,32,29} = 21;


//+
Transfinite Curve{132,103,62,57} = 10 Using Progression 1.15;
Transfinite Curve{129,140,-100,-59,-123,-89,55,86} = 10 Using Progression 1.15;
Transfinite Curve{143,134,112,109} = 10 Using Progression 1.15;


Transfinite Curve{-72,-33,-12,-4} = 10 Using Progression 1.15;
Transfinite Curve{-73,-97,34,-9,53,19,-2,-14} = 10 Using Progression 1.15;
Transfinite Curve{-120,-84,-43,-38} = 10 Using Progression 1.15;


//+
Transfinite Curve{-69,-30,-6,-1,-74,-35,-11,-3,-107,-66,-27,-22,-131,-102,-61,-56} = 10 Using Progression 1.15;


Transfinite Curve{117,81,40,37,121,85,44,39,138,116,80,77,144,135,113,110} = 10 Using Progression 1.15;


//+
Transfinite Curve{118,95,70,68,119,136,-96,-75,-125,-104,71,106,142,139,128,130} = 10 Using Progression 1.15;


Transfinite Curve{-41,-17,-7,-6,-42,-78,18,-8,-47,-24,-5,-10,-26,-111,-88,-58,-60} = 10 Using Progression 1.15;
//+
Transfinite Surface "*";
//+
Recombine Surface "*";
Transfinite Volume "*";
