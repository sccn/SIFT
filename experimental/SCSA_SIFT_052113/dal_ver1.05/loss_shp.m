function [floss, gloss]=loss_shp(zz, bb)

floss = -sum(log(sech(bb-zz)./pi));
gloss = tanh(zz-bb);

% flossp = floss

% keyboard