GFORTRAN module version '6' created from materials.f90 on Mon Aug 13 13:55:06 2012
MD5:d6a9c3cd2f0b34f6add5ad77bb0526a8 -- If you edit this, you'll get what you deserve.

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
(UNKNOWN 0 0 0 UNKNOWN ()) 9 0 (10 11 12) () 0 () () () 0 0)
13 'load_isotope' 'materials' 'load_isotope' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 14 0 (15 16 17 18 19 20 21) () 0 () () () 0 0)
22 'load_source' 'materials' 'load_source' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 23 0
(24 25 26) () 0 () () () 0 0)
27 'material_type' 'materials' 'material_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0
0 UNKNOWN ()) 0 0 () () 0 ((28 'source' (DERIVED 29 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 29 0 0 DERIVED ()) 0 (((NULL (REAL 8
0 0 REAL ()) 0) ()) (() ())) ())) (30 'isotopes' (DERIVED 31 0 0 DERIVED
()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS (STRUCTURE (DERIVED 31
0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (REAL 8 0 0
REAL ()) 0) ()) ((NULL (REAL 8 0 0 REAL ()) 0) ()) ((NULL (REAL 8 0 0
REAL ()) 0) ()) (() ()) (() ()) ((STRUCTURE (DERIVED 32 0 0 DERIVED ())
0 ((() ()) (() ()) ((NULL (REAL 8 0 0 REAL ()) 0) ()) ((NULL (REAL 8 0 0
REAL ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ()) (() ()) ((NULL
(REAL 8 0 0 REAL ()) 0) ()) ((NULL (REAL 8 0 0 REAL ()) 0) ()) ((NULL (
REAL 8 0 0 REAL ()) 0) ()) (() ()) (() ()) (() ()) (() ())) ())) (33
'nisotopes' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (34 'curr_iso' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (35 'npts' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (36 'e_width' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (37
'e_min' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (38 'e_max' (REAL 8 0 0 REAL ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (39 'vol' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (40
'totalxs' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (41 'scattxs' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () ()
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (42 'absorxs' (REAL 8 0 0 REAL
()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (43
'captuxs' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (44 'fissixs' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () ()
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (45 'transxs' (REAL 8 0 0 REAL
()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (46
'xs_total_brdn' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (47 'xs_capt_brdn' (REAL 8 0 0 REAL ()) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (48 'xs_scat_brdn' (REAL 8 0 0
REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (49
'xs_fiss_brdn' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (50 'xs_absb_brdn' (REAL 8 0 0 REAL ()) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (51 'xs_tran_brdn' (REAL 8 0 0
REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ())) PUBLIC (
() () () ()) () 0 0 61897535)
52 'setup_material' 'materials' 'setup_material' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 53 0 (54 55 56 57 58) () 0 () () () 0 0)
29 'source_type' 'materials' 'source_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((59 'e' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE
DIMENSION) UNKNOWN-ACCESS ()) (60 'cdf_width' (REAL 8 0 0 REAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 85089771)
32 'thermal_type' 'materials' 'thermal_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((61 'ktsize' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (62
'cdfsize' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (63 'ktvec' (REAL 8
0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ())
(64 'erat' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (65 'cdf_width' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()))
PUBLIC (() () () ()) () 0 0 35757205)
31 'iso_type' 'materials' 'iso_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOC_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((66 'n' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (67 'a' (REAL 8 0 0
REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (68 'alpha' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (69
'mubar' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (70 'xs_capt' (REAL 8 0 0 REAL ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (71 'xs_scat' (
REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ())
(72 'xs_fiss' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (73 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '255'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (74 'thermal' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (75 'thermal_lib' (DERIVED 32 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 32 0 0 DERIVED ()) 0 ((() ()) (
() ()) ((NULL (REAL 8 0 0 REAL ()) 0) ()) ((NULL (REAL 8 0 0 REAL ()) 0)
()) (() ())) ())) (76 'alpha_mb' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (77
'xs_capt_brdn' (REAL 8 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (78 'xs_scat_brdn'
(REAL 8 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (79 'xs_fiss_brdn' (REAL 8 0 0 REAL ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (80 'engy_capt' (REAL 8 0 0 REAL ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (81 'engy_scat' (REAL 8 0 0
REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION) UNKNOWN-ACCESS ()) (82
'engy_fiss' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 ALLOCATABLE DIMENSION)
UNKNOWN-ACCESS ()) (83 'capt_size' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (84 'scat_size' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (85 'fiss_size' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (86 'doppler' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 1686055)
54 'this' '' 'this' 53 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 27 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
55 'emin' '' 'emin' 53 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
56 'emax' '' 'emax' 53 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
57 'nisotopes' '' 'nisotopes' 53 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
58 'vol' '' 'vol' 53 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
15 'this' '' 'this' 14 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 TARGET DUMMY) (DERIVED 27 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
16 'n' '' 'n' 14 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
17 'a' '' 'a' 14 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
18 'path' '' 'path' 14 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '255'))) 0 0 () () 0 () () () 0 0)
19 'thermal' '' 'thermal' 14 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
20 'name' '' 'name' 14 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '255'))) 0 0 () () 0 () () () 0 0)
21 'doppler' '' 'doppler' 14 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
11 'e' '' 'e' 9 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
10 'this' '' 'this' 9 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 TARGET DUMMY) (DERIVED 27 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
12 'sample_per_xs' '' 'sample_per_xs' 9 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
() 0 () () () 0 0)
4 'this' '' 'this' 3 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 TARGET DUMMY) (DERIVED 27 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
7 'this' '' 'this' 6 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 27 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
24 'this' '' 'this' 23 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 27 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
25 'source_type' '' 'source_type' 23 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
() 0 () () () 0 0)
26 'source_path' '' 'source_path' 23 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '255'))) 0 0 () () 0 () () () 0 0)
)

('compute_macroxs' 0 2 'deallocate_material' 0 5 'doppler_broaden_xs' 0
8 'load_isotope' 0 13 'load_source' 0 22 'material_type' 0 27
'setup_material' 0 52)
