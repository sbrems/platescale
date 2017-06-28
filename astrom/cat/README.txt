The catalog is not formally published yet, so please don't pass it on to other people, or use it for 47Tuc science papers.
It is OK if you or your collaborators want to use it for astrometric-calibration purposes.

For high-precision calibrations, I suggest you use the corrected coordinates in the astrometrically-flat pixel-based frame.
These can be obtained as:

X_cor=x_M(col. 8)-Delta_x(col. 10)
Y_cor=y_M(col. 9)-Delta_y(col. 11)

where x_M,y_M are based on the non-CTE-corrected positions of the original Ata Sarajedini's catalogs, and the Delta_x,Delta_y corrections are derived through iteration in the PM-derivation process.

The reference time of my master frame is 2006.197646 (year, =53808.14079 mjd). The pixel scale of my master frame is 40 mas/pixel.
