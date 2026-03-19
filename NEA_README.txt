Data obtained from the Lowell Observatory website
https://asteroid.lowell.edu/astorb/
Then pre-selected and reformatted using nea_pre_processing.py

File Structure from the website:

Here are two sample records (with one line of parameter numbers above and three lines of column counts below):

(1)    (2)                (3)             (4)    (5)  (6)  (7)   (8)
     1 Ceres              E. Bowell        3.34  0.12 0.72 913.0 G?
  1693 Hertzsprung        E. Bowell       10.97  0.15 0.74  39.5 C
         0         0         0         0         0         0         0
         1         2         3         4         5         6         7
1234567890123456789012345678901234567890123456789012345678901234567890

   (9)                   (10) (11)  (12)     (13)       (14)       (15)
   0   0   0   0   0   0 56959 4750 19960427  80.477333  71.802404  80
   0   0   0   0   0   0 20972   25 19960427 322.276332 234.698906  70
         0         0         1         1         1         1         1
         8         9         0         1         2         3         4
1234567890123456789012345678901234567890123456789012345678901234567890

        (16)      (17)       (18)        (19)     (20)     (21)    (22)
.659857 10.600303 0.07604100   2.76788714 19960414 2.3E-02  1.4E-04 19
.393559 11.942428 0.27460300   2.79629204 19950513 9.0E-01  7.9E-03 19
         1         1         1         1         1         2         2
         5         6         7         8         9         0         1
1234567890123456789012345678901234567890123456789012345678901234567890

      (23)             (24)             (25)
960416 2.7E-02 19960530 3.1E-02 20040111 3.1E-02 20040111
960416 1.2E+00 19960610 1.3E+00 20010812 9.0E-01 20010813
         2         2         2         2         2
         2         3         4         5         6
123456789012345678901234567890123456789012345678901234567
                    

A FORTRAN format statement for reading a record in astorb.dat is (revised 2018/01/02 to allow for an inclination over 99.99999 deg):

A6,1X,A18,1X,A15,1X,A5,1X,F5.2,1X,A4,1X,A5,1X,A4,1X,6I4,1X,
2I5,1X,I4,2I2.2,1X,2(F10.6,1X),F10.6,F10.6,1X,F10.8,F13.8,1X,I4,2I2.2,1X,F7.2,1X,F8.2,1X,I4,2I2,3(1X,F7.2,1X,I4,2I2)

Note that some numerical data (e.g., asteroid number) are encoded as character variables. You may need to decode them. Also, to allow for semi-major axes greater than 1000 AU, a change was made (in 2018) to the output for semi-major axis. If the semi-major axis is less than 1000, the value is written as F13.8; if equal to or larger than 1000, the value is written as F13.7. In either case, the value can be read as F13.8. Note that with these two changes, there may or may not be a space preceeding the inclination and/or the semi-major axis.

Parameters are:
Parameter 	Description
(1) 	Asteroid number (blank if unnumbered).
(2) 	Name or preliminary designation.
(3) 	Orbit computer.
(4) 	Absolute magnitude H, mag [see E. Bowell et al., pp. 549-554, in "Asteroids II", R. P. Binzel et al. (eds.), The University of Arizona Press, Tucson, 1989 and more recent Minor Planet Circulars]. Note that H may be given to 2 decimal places (e.g., 13.41), 1 decimal place (13.4) or as an integer (13), depending on its estimated accuracy. H is given to two decimal places for all unnumbered asteroids, even though it may be very poorly known.
(5) 	Slope parameter G ( ibid.).
(6) 	Color index B-V, mag (blank if unknown; see E. F. Tedesco, pp. 1090-1138, op. cit. ).
(7) 	IRAS diameter, km (blank if unknown; see E. F. Tedesco et al., pp. 1151-1161, op.cit.).
(8) 	IRAS Taxonomic classification (blank if unknown; ibid.).
(9) 	Six integer codes (see table of explanation below). Note that not all codes have been correctly computed.
(10) 	Orbital arc, days, spanned by observations used in orbit computation.
(11) 	Number of observations used in orbit computation.
(12) 	Epoch of osculation, yyyymmdd (TDT). The epoch is the Julian date ending in 00.5 nearest the date the file was created. Thus, as the file is updated, epochs will succeed each other at 100-day intervals on or after Julian dates ending in 50.5 (19980328, 19980706, 19981014, 19990122,...)
(13) 	Mean anomaly, deg.
(14) 	Argument of perihelion, deg (J2000.0).
(15) 	Longitude of ascending node, deg (J2000.0).
(16) 	Inclination, deg (J2000.0).
(17) 	Eccentricity.
(18) 	Semimajor axis, AU.
(19) 	Date of orbit computation, yymmdd (MST, = UTC - 7 hr).
(20) 	Absolute value of the current 1-σ ephemeris uncertainty (CEU), arcsec.
(21) 	Rate of change of CEU, arcsec/day.
(22) 	Date of CEU, yyyymmdd (0 hr UT).
(23) 	Next peak ephemeris uncertainty (PEU), arcsec, from date of CEU, and date of its occurrence, yyyymmdd.
(24) 	Greatest PEU, arcsec, in 10 years from date of CEU, and date of its occurrence, yyyymmdd.
(25) 	Greatest PEU, arcsec, in 10 years from date of next PEU, and date of its occurrence, yyyymmdd, if two observations (of accuracy equal to that of the observations currently included in the orbit--typically ± 1 arcsec) were to be made on the date of the next PEU [parameter (23)].

The meanings of the six integer codes [parameter (9)] are as follows (reference to "type 6:7", for example, means code 6, value 7):
Code 	Value 	Explanation
1 		Planet-crossing asteroids.
		Note: Because some orbits are very poor (or erroneously linked), there may be errors in assignment of these parameter values.
	1 	Aten asteroids (a < 1.0AU).
	2 	Apollo asteroids (a > 1.0AU 0 < q < 1.0).
	4 	Amor asteroids (a > 1.0167AU; 1.0167 < q <
	8 	Mars crossers (1.3 < q < 1.6660AU).
	16 	Outer-planet crossers (excluding Jupiter and Mars Trojans). Asteroids that cross or pass into the heliocentric distance zones between the perihelion and aphelion distances of Jupiter (4.950 to 5.455 AU), Saturn (9.009 to 10.069 AU), Uranus (18.274 to 20 089 AU), and/or Neptune (29.800 to 30.317 AU).
	n 	Asteroids (excluding Mars and Jupiter Trojans) that cross both inner- and outer-planet orbits. For example, an asteroid having n = 24 crosse the orbits of both Mars (q < 1.6660 AU) and Jupiter (Q > 4.950 AU).
2 		Orbit computation.
	1 	Orbits derived from uncertainly, perhaps erroneously linked observations.
	2 	Eccentricity assumed.
	4 	Eccentricity and semimajor axis assumed.
	8 	Mainly for numbered asteroids, omitted observations have resulted in degradation of a so-called orbit-quality parameter (OQP, see K. Muinonen and E. Bowell, Icarus 104, 255-279, 1993), generally by more than 0.1. The corresponding ephemeris uncertainty has increased by about 25% or more.
	16 	OQP degrades by more than 0.1 if unsubstantiated observations (e.g., one-night apparitions) are omitted.
	32 	Orbit derived from data containing observations not in Minor Planet Center files.
	64 	H is unknown. H = 14 mag assumed.
	128 	Asteroid sought, but not found.
	n 	Sum of preceding entries. For example, n = 3 pertains to an uncertainly linked orbit for which the eccentricity was assumed.
3 		Asteroids observed during the course of major surveys. Our definition includes asteroids that were observed but not discovered during the course of a survey.
	1 	Palomar-Leiden survey (PLS) asteroids.
	2 	Palomar-Leiden T-2 survey asteroids.
	4 	Palomar-Leiden T-3 survey asteroids.
	8 	U.K. Schmidt Telescope-Caltech asteroid survey (UCAS) asteroids.
	16 	Palomar-Leiden T-1 survey asteroids.
	n 	Asteroids observed in more than one survey. For example, n = 3 denotes an asteroid observed in both the PLS and T-2 surveys.
4 		Minor Planet Center (MPC) critical-list numbered asteroids.
	1 	Lost asteroid.
	2 	Asteroids observed at only two apparitions.
	3 	Asteroids observed at only three apparitions.
	4 	Asteroids observed at four or more apparitions, last more than ten years ago.
	5 	Asteroids observed at four or more apparitions, only one night in last ten years.
	6 	Other poorly observed asteroids observed at four or more apparitions.
	7 	Absolute magnitude poorly known (not on MPC critical-list).
5 		Lowell Observatory and related discoveries
	1 	Asteroids discovered by E. Bowell.
	2 	Non-Bowell discoveries from Lowell search programs.
	3 	Sum of preceding entries. n = 3 pertains to an asteroid discovered jointly by E. Bowell and another person connected with Lowell search programs.
6 		Rank, in decreasing importance, for our collaborative program of astrometry using the transit circle of the U.S. Naval Observatory Flagstaff Station.
	10 	Exceptionally important, to be observed frequently. Principally space mission targets and occultation candidates.
	9 	Asteroids useful for mass determination.
	8 	Asteroids for which one or two additional nights' observation are required to satisfy orbit-update requirements. Asteroids of type 6:7 whose ephemeris uncertainties are between 2 and 5 arcsec within the next ten years or so.
	7 	Bowell unnumbered discoveries whose ephemeris uncertainties are less than 2 arcsec within the next ten years or so. MPC critical-list asteroids.
	6 	Planet-crossers of type 6:5.
	5 	Numbered asteroids whose ephemeris uncertainties are between 2 and 5 arcsec within the next ten years or so. Unnumbered asteroids that should be numberable after one or two more nights' observation.

Note that the codes have not been carefully checked. There are doubtless many errors. 
