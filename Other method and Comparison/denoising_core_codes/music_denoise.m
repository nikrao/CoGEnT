function [x_music,fs_music,cs_music] =  music_denoise(observed, k)
n = length(observed);
fs_music = sort(mod(rootmusic(observed,k)/2/pi,1));
music_basis = exp( 1i*(0:(n-1))'*2*pi*reshape(fs_music,1,k) );
cs_music = music_basis\observed;
x_music = music_basis*cs_music;
end
