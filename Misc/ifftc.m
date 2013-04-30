function [ dataArray ] = ifftc( dataArray )
%IFFTC inverse FFT for centered data 
%
%   Syntax
%
%   A = IFFTS( F )
%
%   Returns IFFT of array F after applying the correct ifftshifts
%

dataArray = fftshift( ifftn( ( dataArray ) ) ) ; 

end