import os
from textwrap import dedent as ddent
import shutil
from astropy.time import Time

def eq1file(orbit, tdir):
    '''
    Parameters
    ----------
    el_jpl : `~numpy.ndarray` (N, 18)
        Orbital elements as determined by JPL Horizons.
    tdir : `tf.TemporaryDirectory`
        Temporary directory object, temp directory for file generation
        
    Returns
    ------- 
    temp.eq1: temporary text file to run OrbFit
    '''  
    
    os.makedirs(f'{tdir}/epoch',exist_ok=True)
    
    with open(f"{tdir}/epoch/temp.eq1", "w") as fp:
        fp.write(ddent(f'''\
        format  = 'OEF2.0'       ! file format
        rectype = 'ML'           ! record type (1L/ML)
        refsys  = ECLM J2000     ! default reference system
        END_OF_HEADER
        temp
        ! Cometary elements: q, e, i, long. node, arg. peric., pericenter time\n'''))
        car = orbit.cartesian
        st1=" CAR   %.16E  %.16E   %.16E  %.16E  %.16E   %.16E\n" % (car[0][0], car[0][1], car[0][2], car[0][3], car[0][4], car[0][5])
        #fp.write(" COM   %.15E  %.15f   %.15f  %.15f  %.15f   %.15f\n" % (el_jpl['q'][0], el_jpl['e'][0], el_jpl['incl'][0], el_jpl['Omega'][0], el_jpl['w'][0], ap.time.Time(el_jpl['Tp_jd'][0], format='jd').mjd))
        fp.write(st1)
        fp.write(" MJD     %.9f TDT\n" % (orbit.epochs[0].tt.mjd)) # tdt
        fp.write(" MAG  %.3f  %.3f\n" % (orbit.H[0], orbit.G[0]))
        fp.write(ddent(f'''\
        ! Non-grav parameters: model used, actual number in use, dimension
         LSP   0  0    6
        ! RMS    1.68232E-08   5.88900E-08   8.22688E-08   5.34469E-08   7.56890E-08   7.13398E-06
        ! EIG   1.07579E-08   5.31009E-08   5.91864E-08   7.57356E-08   7.86169E-08   1.27486E-07
        ! WEA   0.08651  -0.02289  -0.24412  -0.03357  -0.02078  -0.96480
         COV   2.830210138102556E-16  3.543122024213312E-16 -2.603292682702056E-16
         COV  -4.992042214900484E-18 -1.052690180196314E-18 -7.873861865190710E-14
         COV   3.468027286557968E-15  6.878077752471183E-17  1.886511729787680E-17
         COV   6.689670038864485E-17  1.808279351538482E-14  6.768149177265766E-15
         COV   7.159243161286040E-17  1.248926483233068E-16  1.357728291186093E-13
         COV   2.856568560436748E-15  2.588049637167598E-16  2.529981071526617E-14
         COV   5.728829671236270E-15  1.056596023015451E-14  5.089368128905899E-11
         NOR   8.462990106959648E+15 -9.345934921051774E+14  6.961302078833404E+13
         NOR  -9.766026206616650E+13 -9.148695418092123E+12  1.329006055003970E+13
         NOR   3.921580648140324E+14 -8.566206317904612E+12  1.006265833999790E+13
         NOR  -2.128841368531985E+12 -1.566971456817283E+12  1.567202493495569E+14
         NOR  -7.910041612493922E+11 -2.702958007388599E+12 -3.063965034542373E+11
         NOR   3.541407046591562E+14 -1.551670529664669E+13 -3.253830675316872E+11
         NOR   1.754031538722264E+14 -3.488851624201696E+10  4.175326401599722E+10
         '''))
    
    return

def eq1filePROP(orbit, tdir):
    '''
    Parameters
    ----------
    el_jpl : `~numpy.ndarray` (N, 18)
        Orbital elements as determined by JPL Horizons.
    tdir : `tf.TemporaryDirectory`
        Temporary directory object, temp directory for file generation
        
    Returns
    ------- 
    temp.eq1: temporary text file to run OrbFit
    '''  
    
    os.makedirs(f'{tdir}/epoch',exist_ok=True)
    
    with open(f"{tdir}/epoch/temp.eq1", "w") as fp:
        fp.write(ddent(f'''\
        format  = 'OEF2.0'       ! file format
        rectype = 'ML'           ! record type (1L/ML)
        refsys  = ECLM J2000     ! default reference system
        END_OF_HEADER
        temp
        ! Cometary elements: q, e, i, long. node, arg. peric., pericenter time\n'''))
        car = orbit
        st1=" CAR   %.16E  %.16E   %.16E  %.16E  %.16E   %.16E\n" % (car['x'][0], car['y'][1], car['z'][2], car['vx'][3], car['vy'][4], car['vz'][5])
        #fp.write(" COM   %.15E  %.15f   %.15f  %.15f  %.15f   %.15f\n" % (el_jpl['q'][0], el_jpl['e'][0], el_jpl['incl'][0], el_jpl['Omega'][0], el_jpl['w'][0], ap.time.Time(el_jpl['Tp_jd'][0], format='jd').mjd))
        fp.write(st1)
        fp.write(" MJD     %.9f TDT\n" % (orbit.epochs[0].tt.mjd))
        fp.write(" MAG  %.3f  %.3f\n" % (orbit.H[0], orbit.G[0]))
        fp.write(ddent(f'''\
        ! Non-grav parameters: model used, actual number in use, dimension
         LSP   0  0    6
        ! RMS    1.68232E-08   5.88900E-08   8.22688E-08   5.34469E-08   7.56890E-08   7.13398E-06
        ! EIG   1.07579E-08   5.31009E-08   5.91864E-08   7.57356E-08   7.86169E-08   1.27486E-07
        ! WEA   0.08651  -0.02289  -0.24412  -0.03357  -0.02078  -0.96480'''))
        #  COV   2.830210138102556E-16  3.543122024213312E-16 -2.603292682702056E-16
        #  COV  -4.992042214900484E-18 -1.052690180196314E-18 -7.873861865190710E-14
        #  COV   3.468027286557968E-15  6.878077752471183E-17  1.886511729787680E-17
        #  COV   6.689670038864485E-17  1.808279351538482E-14  6.768149177265766E-15
        #  COV   7.159243161286040E-17  1.248926483233068E-16  1.357728291186093E-13
        #  COV   2.856568560436748E-15  2.588049637167598E-16  2.529981071526617E-14
        #  COV   5.728829671236270E-15  1.056596023015451E-14  5.089368128905899E-11
        #  NOR   8.462990106959648E+15 -9.345934921051774E+14  6.961302078833404E+13
        #  NOR  -9.766026206616650E+13 -9.148695418092123E+12  1.329006055003970E+13
        #  NOR   3.921580648140324E+14 -8.566206317904612E+12  1.006265833999790E+13
        #  NOR  -2.128841368531985E+12 -1.566971456817283E+12  1.567202493495569E+14
        #  NOR  -7.910041612493922E+11 -2.702958007388599E+12 -3.063965034542373E+11
        #  NOR   3.541407046591562E+14 -1.551670529664669E+13 -3.253830675316872E+11
        #  NOR   1.754031538722264E+14 -3.488851624201696E+10  4.175326401599722E+10
        #  '''))
    
    return

# Writes temporary .fop file for OrbFit
def fopfile(tdir,obsc=''):
    '''
    Parameters
    ----------
    tdir : `tf.TemporaryDirectory`
        Temporary directory object, temp directory for file generation
    
    Returns
    ------- 
    temp.fop: temporary text file to run OrbFit
    '''
        
    with open(f"{tdir}/temp.fop", "w") as fp:
        fp.write(ddent(f'''\
        ! input file for fitobs
        fitobs.
        ! first arc        .astna0='temp'           ! full name, first arc
                .obsdir0='mpcobs/'         ! directory of observs. file, first arc
                .elefi0='epoch/temp.eq1' ! first arc elements file

        ! second arc
        !        .astnap=''            ! full name, second arc
        !        .obsdirp='mpcobs'     ! directory of observs. file, second arc
        ! bizarre  control;
                .ecclim=     1.9999d0    ! max eccentricity for non bizarre orbit
                .samin=      0.3d0       ! min a for non bizarre orbit
                .samax=      2000.d0     ! max a for non bizarre orbit
                .phmin=      0.001d0     ! min q for non bizarre orbit
                .ahmax=     4000.d0      ! max Q for non bizarre orbit
                .error_model='fcct14'     ! error model
        ! ephem.
                .timescale='UTC'
                .obscode={obsc}
        propag.
                .iast=17            ! 0=no asteroids with mass, n=no. of massive asteroids (def=0)
                .filbe='AST17'      ! name of the asteroid ephemerides file (def='CPV')
                .npoint=600         ! minimum number of data points for a deep close appr (def=100)
                .dmea=0.2d0         ! min. distance for control of close-app. to the Earth only (def=0.1)
                .dter=0.05d0        ! min. distance for control of close-app.
                                    ! to terrestrial planets (MVM)(def=0.1)
                .yark_exp=2.d0      ! A2/r^yark_exp model (def=2)
                .ngr_opt=.FALSE.     ! read options for non-gravitational perturbations from the option file 
                .irel=2             ! 0=newtonian 1=gen. relativity, sun 2=gen. rel. all planets
                                    !          (def=0, 1 for NEA, 2 for radar)
                .iaber=2            ! aberration 0=no 1=yes 2=(def=1)
                .ilun=1             ! 0=no moon 1= yes (def=0, 1 for NEA)
                .iyark=3            ! 0=no Yarkovsky, 1=Yark diurnal, 2=Yark seasonal
                                    !    3=secular nongravitational perturbations (including Yark) (def=0)
                .ipa2m=0           ! 0=no drpa2m, 1=yes spherical direct radiation pressure (def=0)
                .det_drp=2          ! how many parameters to solve: 0=none 1=drpa2m 2=dadt 3=both (def=0)
                .det_outgas=0       ! det outgassing for comets
        

        difcor.

        IERS.
                .extrapolation=.T. ! extrapolation of Earth rotation

        reject.
                .rejopp=.false.    ! reject entire opposition
        '''))
    home = os.environ["HOME"]
    #home = os.path.join(os.path.dirname(os.path.dirname(__file__)), "orbf") 
    shutil.copyfile(home+"/orbfit/lib/AST17.bai", f"{tdir}/AST17.bai")
    shutil.copyfile(home+"/orbfit/lib/AST17.bep", f"{tdir}/AST17.bep")

    return


def propHelp(stop, tdir):
    '''
    Parameters
    ----------
    start : `str`
        Beginning time of integration, UTC.  
    stop : `str`
        End time of integration, UTC.   
    obs : `str`
        Observatory code.
    tdir : `tf.TemporaryDirectory`
        Temporary directory object, temp directory for file generation
        
    Returns
    ------- 
    ast.inp: text file required to bypass OrbFit interactive menu
    '''
    
    with open(f"{tdir}/ast.inp", "w") as fp:
        fp.write(ddent(f'''\
        temp
        5
        1
        1
        {stop.tt.mjd}
        0
        '''))
        
    return

def ephHelp(start, stop, step, obs, tdir):
    '''
    Parameters
    ----------
    start : `str`
        Beginning time of integration, UTC.  
    stop : `str`
        End time of integration, UTC.   
    obs : `str`
        Observatory code.
    tdir : `tf.TemporaryDirectory`
        Temporary directory object, temp directory for file generation
        
    Returns
    ------- 
    ast.inp: text file required to bypass OrbFit interactive menu
    '''
    
    with open(f"{tdir}/ast.inp", "w") as fp:
        fp.write(ddent(f'''\
        temp
        6
        6
        {start.tt.mjd}
        {stop.tt.mjd}
        {step:0.0f}
        {obs}
        0
        '''))
        
    return

def fopOD(tdir):
    '''
    Parameters
    ----------
    tdir : `tf.TemporaryDirectory`
        Temporary directory object, temp directory for file generation
    
    Returns
    ------- 
    temp.fop: temporary text file to run OrbFit
    '''
        
    with open(f"{tdir}/temp.fop", "w") as fp:
        fp.write(ddent('''\
        ! input file for fitobs
        fitobs.
        ! first arc        .astna0='temp'           ! full name, first arc
                .obsdir0='obs/'         ! directory of observs. file, first arc
                .elefi0='epoch/temp.eq1' ! first arc elements file

        ! second arc
        !        .astnap=''            ! full name, second arc
        !        .obsdirp='mpcobs'     ! directory of observs. file, second arc
        ! bizarre  control;
                .ecclim=     1.9999d0    ! max eccentricity for non bizarre orbit
                .samin=      0.3d0       ! min a for non bizarre orbit
                .samax=      2000.d0     ! max a for non bizarre orbit
                .phmin=      0.001d0     ! min q for non bizarre orbit
                .ahmax=     4000.d0      ! max Q for non bizarre orbit
                .error_model='fcct14'     ! error model
        ! ephem.
                .timescale='UTC'
        propag.
                .iast=17            ! 0=no asteroids with mass, n=no. of massive asteroids (def=0)
                .filbe='AST17'      ! name of the asteroid ephemerides file (def='CPV')
                .npoint=600         ! minimum number of data points for a deep close appr (def=100)
                .dmea=0.2d0         ! min. distance for control of close-app. to the Earth only (def=0.1)
                .dter=0.05d0        ! min. distance for control of close-app.
                                    ! to terrestrial planets (MVM)(def=0.1)
                .yark_exp=2.d0      ! A2/r^yark_exp model (def=2)
                .ngr_opt=.FALSE.     ! read options for non-gravitational perturbations from the option file
                .irel=2             ! 0=newtonian 1=gen. relativity, sun 2=gen. rel. all planets
                                    !          (def=0, 1 for NEA, 2 for radar)
                .iaber=2            ! aberration 0=no 1=yes 2=(def=1)
                .ilun=1             ! 0=no moon 1= yes (def=0, 1 for NEA)
                .iyark=3            ! 0=no Yarkovsky, 1=Yark diurnal, 2=Yark seasonal
                                    !    3=secular nongravitational perturbations (including Yark) (def=0)
                .ipa2m=0           ! 0=no drpa2m, 1=yes spherical direct radiation pressure (def=0)
                .det_drp=2          ! how many parameters to solve: 0=none 1=drpa2m 2=dadt 3=both (def=0)
                .det_outgas=0       ! det outgassing for comets

        difcor.

        IERS.
                .extrapolation=.T. ! extrapolation of Earth rotation

        reject.
                .rejopp=.false.    ! reject entire opposition
        '''))
    home = os.environ["HOME"]
    #home = os.path.join(os.path.dirname(os.path.dirname(__file__)), "orbf") 
    shutil.copyfile(home+"/orbfit/lib/AST17.bai", f"{tdir}/AST17.bai")
    shutil.copyfile(home+"/orbfit/lib/AST17.bep", f"{tdir}/AST17.bep")

    return

def fopOD2(tdir):
    '''
    Parameters
    ----------
    tdir : `tf.TemporaryDirectory`
        Temporary directory object, temp directory for file generation
    
    Returns
    ------- 
    temp.fop: temporary text file to run OrbFit
    '''
        
    with open(f"{tdir}/temp.fop", "w") as fp:
        fp.write(ddent('''\
        ! input file for fitobs
        fitobs.
        ! first arc        .astna0='temp'           ! full name, first arc
                .obsdir0='obs1/'         ! directory of observs. file, first arc
                .elefi0='epoch/temp.eq1' ! first arc elements file

        ! second arc
        !        .astnap=''            ! full name, second arc
        !        .obsdirp='mpcobs'     ! directory of observs. file, second arc
        ! bizarre  control;
                .ecclim=     1.9999d0    ! max eccentricity for non bizarre orbit
                .samin=      0.3d0       ! min a for non bizarre orbit
                .samax=      2000.d0     ! max a for non bizarre orbit
                .phmin=      0.001d0     ! min q for non bizarre orbit
                .ahmax=     4000.d0      ! max Q for non bizarre orbit
                .error_model='fcct14'     ! error model
        ! ephem.
                .timescale='UTC'
        propag.
                .iast=17            ! 0=no asteroids with mass, n=no. of massive asteroids (def=0)
                .filbe='AST17'      ! name of the asteroid ephemerides file (def='CPV')
                .npoint=600         ! minimum number of data points for a deep close appr (def=100)
                .dmea=0.2d0         ! min. distance for control of close-app. to the Earth only (def=0.1)
                .dter=0.05d0        ! min. distance for control of close-app.
                                    ! to terrestrial planets (MVM)(def=0.1)
                .yark_exp=2.d0      ! A2/r^yark_exp model (def=2)
                .ngr_opt=.FALSE.     ! read options for non-gravitational perturbations from the option file
                .irel=2             ! 0=newtonian 1=gen. relativity, sun 2=gen. rel. all planets
                                    !          (def=0, 1 for NEA, 2 for radar)
                .iaber=2            ! aberration 0=no 1=yes 2=(def=1)
                .ilun=1             ! 0=no moon 1= yes (def=0, 1 for NEA)
                .iyark=3            ! 0=no Yarkovsky, 1=Yark diurnal, 2=Yark seasonal
                                    !    3=secular nongravitational perturbations (including Yark) (def=0)
                .ipa2m=0           ! 0=no drpa2m, 1=yes spherical direct radiation pressure (def=0)
                .det_drp=2          ! how many parameters to solve: 0=none 1=drpa2m 2=dadt 3=both (def=0)
                .det_outgas=0       ! det outgassing for comets

        difcor.

        IERS.
                .extrapolation=.T. ! extrapolation of Earth rotation

        reject.
                .rejopp=.false.    ! reject entire opposition
        '''))
    home = os.environ["HOME"]
    #home = os.path.join(os.path.dirname(os.path.dirname(__file__)), "orbf") 
    shutil.copyfile(home+"/orbfit/lib/AST17.bai", f"{tdir}/AST17.bai")
    shutil.copyfile(home+"/orbfit/lib/AST17.bep", f"{tdir}/AST17.bep")

    return

#FOR OD PART 2
def eqOD(orbit, tdir):
    '''
    Parameters
    ----------
    el_jpl : `~numpy.ndarray` (N, 18)
        Orbital elements as determined by JPL Horizons.
    tdir : `tf.TemporaryDirectory`
        Temporary directory object, temp directory for file generation
        
    Returns
    ------- 
    temp.eq1: temporary text file to run OrbFit
    '''  
    
    os.makedirs(f'{tdir}/epoch',exist_ok=True)
    
    with open(f"{tdir}/epoch/temp.eq1", "w") as fp:
        fp.write(ddent(f'''\
        format  = 'OEF2.0'       ! file format
        rectype = 'ML'           ! record type (1L/ML)
        refsys  = ECLM J2000     ! default reference system
        END_OF_HEADER
        temp
        ! Cometary elements: q, e, i, long. node, arg. peric., pericenter time\n'''))
        car = orbit
        st1=" KEP   %.16E  %.16E   %.16E  %.16E  %.16E   %.16E\n" % (car[0], car[1], car[2], car[3], car[4], car[5])
        #fp.write(" COM   %.15E  %.15f   %.15f  %.15f  %.15f   %.15f\n" % (el_jpl['q'][0], el_jpl['e'][0], el_jpl['incl'][0], el_jpl['Omega'][0], el_jpl['w'][0], ap.time.Time(el_jpl['Tp_jd'][0], format='jd').mjd))
        fp.write(st1)
        fp.write(" MJD     %.9f TDT\n" % (Time(orbit[6],format='mjd',scale='tt').tt.mjd)) # CHECK THIS
        #fp.write(" MAG  %.3f  %.3f\n" % (orbit.H[0], orbit.G[0]))
        fp.write(ddent(f'''\
        ! Non-grav parameters: model used, actual number in use, dimension
         LSP   0  0    6
        ! RMS    1.68232E-08   5.88900E-08   8.22688E-08   5.34469E-08   7.56890E-08   7.13398E-06
        ! EIG   1.07579E-08   5.31009E-08   5.91864E-08   7.57356E-08   7.86169E-08   1.27486E-07
        ! WEA   0.08651  -0.02289  -0.24412  -0.03357  -0.02078  -0.96480
         COV   2.830210138102556E-16  3.543122024213312E-16 -2.603292682702056E-16
         COV  -4.992042214900484E-18 -1.052690180196314E-18 -7.873861865190710E-14
         COV   3.468027286557968E-15  6.878077752471183E-17  1.886511729787680E-17
         COV   6.689670038864485E-17  1.808279351538482E-14  6.768149177265766E-15
         COV   7.159243161286040E-17  1.248926483233068E-16  1.357728291186093E-13
         COV   2.856568560436748E-15  2.588049637167598E-16  2.529981071526617E-14
         COV   5.728829671236270E-15  1.056596023015451E-14  5.089368128905899E-11
         NOR   8.462990106959648E+15 -9.345934921051774E+14  6.961302078833404E+13
         NOR  -9.766026206616650E+13 -9.148695418092123E+12  1.329006055003970E+13
         NOR   3.921580648140324E+14 -8.566206317904612E+12  1.006265833999790E+13
         NOR  -2.128841368531985E+12 -1.566971456817283E+12  1.567202493495569E+14
         NOR  -7.910041612493922E+11 -2.702958007388599E+12 -3.063965034542373E+11
         NOR   3.541407046591562E+14 -1.551670529664669E+13 -3.253830675316872E+11
         NOR   1.754031538722264E+14 -3.488851624201696E+10  4.175326401599722E+10
         '''))
    
    return
# FOR Initial Orbit State guess
def od1Help(tdir):
    '''
    Parameters
    ----------
    start : `str`
        Beginning time of integration, UTC.  
    stop : `str`
        End time of integration, UTC.   
    obs : `str`
        Observatory code.
    tdir : `tf.TemporaryDirectory`
        Temporary directory object, temp directory for file generation
        
    Returns
    ------- 
    ast.inp: text file required to bypass OrbFit interactive menu
    '''
    
    with open(f"{tdir}/ast.inp", "w") as fp:
        fp.write(ddent(f'''\
        temp
        2
        4
        1
        0
        '''))
        
    return
# FOR Least Squares improvements
def od2Help(tdir):
    '''
    Parameters
    ----------
    start : `str`
        Beginning time of integration, UTC.  
    stop : `str`
        End time of integration, UTC.   
    obs : `str`
        Observatory code.
    tdir : `tf.TemporaryDirectory`
        Temporary directory object, temp directory for file generation
        
    Returns
    ------- 
    ast.inp: text file required to bypass OrbFit interactive menu
    '''
    
    with open(f"{tdir}/ast.inp", "w") as fp:
        fp.write(ddent(f'''\
        temp
        8
        3
        1
        y
        3
        1
        0
        '''))
        
    return