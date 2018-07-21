settings.outformat = "pdf";
settings.prc = false;
settings.render = 16;
size(5cm,5cm,IgnoreAspect);
import three;
import graph3;

currentprojection=perspective((3.5,-4.5,1.4),up=Z);

real myopacity=1;
real radius=0.05;
pen atompen=white;

draw(unitcube, surfacepen=white+opacity(0.2));
draw(scale3(radius)*unitsphere, surfacepen=atompen);
draw(shift(0,0.5,0.5)*scale3(radius)*unitsphere, surfacepen=atompen);
draw(shift(0.5,0,0.5)*scale3(radius)*unitsphere, surfacepen=atompen);
draw(shift(0.5,0.5,0)*scale3(radius)*unitsphere, surfacepen=atompen);
