GFORTRAN module version '6' created from materials.f90 on Mon May  7 21:38:01 2012
MD5:7def894f7aa574d1fe8fd6fc6d8c7eb4 -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

()

()

()

()

(2 'compute_macroxs' 'materials' 'compute_macroxs' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 3 0 (4) () 0 () () () 0 0)
5 'deallocate_material' 'materials' 'deallocate_material' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 6 0 (7) () 0 () () () 0 0)
8 'doppler_broaden_xs' 'materials' 'doppler_broaden_xs' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 9 0 (10 11) () 0 () () () 0 0)
12 'load_isotope' 'materials' 'load_isotope' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 13 0 (14 15 16 17 18 19 20) () 0 () () () 0 0)
21 'load_source' 'materials' 'load_source' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 22 0
(23 24 25) () 0 () () () 0 0)
26 'material_type' 'materials' 'material_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0
0 UNKNOWN ()) 0 0 () () 0 ((27 'source' (DERIVED 28 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 28 0 0 DERIVED ()) 0 (((NULL (REAL 8
0 0 REAL ()) 0) ()) (() ())) ())) (29 'isotopes' (DERIVED 30 0 0 DERIVED
()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS (STRUCTURE (DERIVED 30
0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (REAL 8 0 0
REAL ()) 0) ()) ((NULL (REAL 8 0 0 REAL ()) 0) ()) ((NULL (REAL 8 0 0
REAL ()) 0) ()) (() ()) (() ()) ((STRUCTURE (DERIVED 31 0 0 DERIVED ())
0 ((() ()) (() ()) ((NULL (REAL 8 0 0 REAL ()) 0) ()) ((NULL (REAL 8 0 0
REAL ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ()) (() ()) ((NULL
(REAL 8 0 0 REAL ()) 0) ()) ((NULL (REAL 8 0 0 REAL ()) 0) ()) ((NULL (
REAL 8 0 0 REAL ()) 0) ()) (() ()) (() ()) (() ()) (() ())) ())) (32
'nisotopes' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (33 'curr_iso' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (34 'npts' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (35 'e_width' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (36
'e_min' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (37 'e_max' (REAL 8 0 0 REAL ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (38 'vol' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (39
'totalxs' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (40 'scattxs' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () ()
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (41 'absorxs' (REAL 8 0 0 REAL
()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (42
'captuxs' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (43 'fissixs' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () ()
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (44 'transxs' (REAL 8 0 0 REAL
()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (45
'xs_total_brdn' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (46 'xs_capt_brdn' (REAL 8 0 0 REAL ()) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (47 'xs_scat_brdn' (REAL 8 0 0
REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (48
'xs_fiss_brdn' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (49 'xs_absb_brdn' (REAL 8 0 0 REAL ()) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (50 'xs_tran_brdn' (REAL 8 0 0
REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ())) PUBLIC (
() () () ()) () 0 0 61897535)
51 'setup_material' 'materials' 'setup_material' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 52 0 (53 54 55 56 57) () 0 () () () 0 0)
28 'source_type' 'materials' 'source_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((58 'e' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE
DIMENSION) UNKNOWN-ACCESS ()) (59 'cdf_width' (REAL 8 0 0 REAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 85089771)
31 'thermal_type' 'materials' 'thermal_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((60 'ktsize' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (61
'cdfsize' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (62 'ktvec' (REAL 8
0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ())
(63 'erat' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (64 'cdf_width' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()))
PUBLIC (() () () ()) () 0 0 35757205)
30 'iso_type' 'materials' 'iso_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((65 'n' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (66 'a' (REAL 8 0 0
REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (67 'alpha' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (68
'mubar' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (69 'xs_capt' (REAL 8 0 0 REAL ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (70 'xs_scat' (
REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ())
(71 'xs_fiss' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (72 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '255'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (73 'thermal' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (74 'thermal_lib' (DERIVED 31 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 31 0 0 DERIVED ()) 0 ((() ()) (
() ()) ((NULL (REAL 8 0 0 REAL ()) 0) ()) ((NULL (REAL 8 0 0 REAL ()) 0)
()) (() ())) ())) (75 'alpha_mb' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (76
'xs_capt_brdn' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (77 'xs_scat_brdn'
(REAL 8 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (78 'xs_fiss_brdn' (REAL 8 0 0 REAL ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (79 'engy_capt' (REAL 8 0 0 REAL ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (80 'engy_scat' (REAL 8 0 0
REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (81
'engy_fiss' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (82 'capt_size' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (83 'scat_size' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (84 'fiss_size' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (85 'doppler' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 1686055)
53 'this' '' 'this' 52 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 26 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
54 'emin' '' 'emin' 52 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
55 'emax' '' 'emax' 52 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
56 'nisotopes' '' 'nisotopes' 52 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
57 'vol' '' 'vol' 52 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
14 'this' '' 'this' 13 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 TARGET DUMMY) (DERIVED 26 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
15 'n' '' 'n' 13 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
16 'a' '' 'a' 13 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
17 'path' '' 'path' 13 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '255'))) 0 0 () () 0 () () () 0 0)
18 'thermal' '' 'thermal' 13 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
19 'name' '' 'name' 13 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '255'))) 0 0 () () 0 () () () 0 0)
20 'doppler' '' 'doppler' 13 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
23 'this' '' 'this' 22 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 26 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
24 'source_type' '' 'source_type' 22 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
() 0 () () () 0 0)
25 'source_path' '' 'source_path' 22 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '255'))) 0 0 () () 0 () () () 0 0)
11 'e' '' 'e' 9 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
10 'this' '' 'this' 9 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 TARGET DUMMY) (DERIVED 26 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
4 'this' '' 'this' 3 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 TARGET DUMMY) (DERIVED 26 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
7 'this' '' 'this' 6 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 26 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
)

('compute_macroxs' 0 2 'deallocate_material' 0 5 'doppler_broaden_xs' 0
8 'load_isotope' 0 12 'load_source' 0 21 'material_type' 0 26
'setup_material' 0 51)
