# this script is for the mid-range SOG modes, here we use the 3DFFT to handle them
# including the following steps
# 1. gridding: interpolate all sources in to uniform real grid
# 2. 3D FFT
# 3. scaling
# 4. 3D iFFT
# 5. gathering: gather all potential from uniform real grid to particles

