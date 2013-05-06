function [ dataArray ] = fftc( dataArray )
%FFTC n-d fft for centered data
%
%   Syntax
%
%   F = FFTC( A )
%
%   Returns FFT of array A after applying the correct fftshifts
%


dataArray = fftn( ifftshift( dataArray ) ) ; 

end