function [LL]=LLDSGE(theta)
global Z

[A,C,D,R]=DSGE_SS(theta);
LL=LL_state_space(A,C,D,R,Z);
