int steps_guess(long long *bs,long long *gs,long long l)
{
  /* l=3: bs=0 gs=0 bench=10898 baseline=10866 */
  /* l=5: bs=0 gs=0 bench=13996 baseline=13874 */
  /* l=7: bs=0 gs=0 bench=19770 baseline=19650 */
  /* l=11: bs=0 gs=0 bench=30392 baseline=30060 */
  /* l=31: bs=0 gs=0 bench=81956 baseline=81678 */
  /* l=43: bs=0 gs=0 bench=111358 baseline=110954 */
  /* l=83: bs=0 gs=0 bench=210558 baseline=210460 */
  /* l=103: bs=0 gs=0 bench=260916 baseline=277756 */
  /* l=107: bs=0 gs=0 bench=271658 baseline=271148 */
  /* l=109: bs=0 gs=0 bench=278950 baseline=275532 */
  /* l=137: bs=0 gs=0 bench=344948 baseline=342378 */
  /* l=199: bs=0 gs=0 bench=499932 baseline=496918 */
  /* l=227: bs=0 gs=0 bench=567720 baseline=566702 */
  if (l <= 227) { *bs = 0; *gs = 0; return 1; }
  /* l=419: bs=14 gs=7 bench=927284 baseline=1039492 */
  if (l <= 419) { *bs = 14; *gs = 7; return 1; }
  /* l=491: bs=14 gs=8 bench=1048384 baseline=1222546 */
  if (l <= 491) { *bs = 14; *gs = 8; return 1; }
  /* l=569: bs=14 gs=10 bench=1164240 baseline=1432688 */
  if (l <= 569) { *bs = 14; *gs = 10; return 1; }
  /* l=631: bs=14 gs=11 bench=1292240 baseline=1569390 */
  if (l <= 631) { *bs = 14; *gs = 11; return 1; }
  /* l=677: bs=16 gs=10 bench=1354698 baseline=1699296 */
  if (l <= 677) { *bs = 16; *gs = 10; return 1; }
  /* l=751: bs=16 gs=11 bench=1491070 baseline=1880982 */
  if (l <= 751) { *bs = 16; *gs = 11; return 1; }
  /* l=827: bs=16 gs=12 bench=1627010 baseline=2079036 */
  if (l <= 827) { *bs = 16; *gs = 12; return 1; }
  /* l=857: bs=16 gs=13 bench=1664550 baseline=2171694 */
  /* l=859: bs=16 gs=13 bench=1666950 baseline=2165976 */
  if (l <= 859) { *bs = 16; *gs = 13; return 1; }
  /* l=883: bs=20 gs=11 bench=1714210 baseline=2206862 */
  if (l <= 883) { *bs = 20; *gs = 11; return 1; }
  /* l=1019: bs=18 gs=14 bench=1914868 baseline=2571866 */
  if (l <= 1019) { *bs = 18; *gs = 14; return 1; }
  /* l=1171: bs=22 gs=13 bench=2141964 baseline=2909050 */
  if (l <= 1171) { *bs = 22; *gs = 13; return 1; }
  /* l=1879: bs=30 gs=15 bench=3003372 baseline=4850172 */
  if (l <= 1879) { *bs = 30; *gs = 15; return 1; }
  /* l=2713: bs=32 gs=21 bench=3999778 baseline=6837360 */
  if (l <= 2713) { *bs = 32; *gs = 21; return 1; }
  /* l=3691: bs=38 gs=24 bench=5223138 baseline=9152296 */
  if (l <= 3691) { *bs = 38; *gs = 24; return 1; }
  /* l=4019: bs=40 gs=25 bench=5592828 baseline=9977994 */
  if (l <= 4019) { *bs = 40; *gs = 25; return 1; }
  /* l=4283: bs=38 gs=28 bench=5876690 baseline=10815114 */
  if (l <= 4283) { *bs = 38; *gs = 28; return 1; }
  /* l=6983: bs=62 gs=28 bench=8551860 baseline=17354004 */
  if (l <= 6983) { *bs = 62; *gs = 28; return 1; }
  return 0;
}
