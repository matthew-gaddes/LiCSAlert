#!/bin/bash

# Note - this will only work if you can login to Jasmin
# which requires running a command such as:
# eval $(ssh-agent -s); ssh-add ~/.ssh/id_rsa_jasmin
# after you have setup a key pair with Jasmin.


rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/alcedo_106A_09048_000909 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/alcedo_107A_09018_000205 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/alcedo_128D_09016_110500 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/cerro_azul_106A_09048_000909 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/cerro_azul_106A_09090_000404 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/cerro_azul_128D_09016_110500 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/puyehue_cordon_caulle_083D_13027_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/puyehue_cordon_caulle_164A_12971_071312 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/puyehue_cordon_caulle_164A_13146_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cordon_del_azufre_047A_11581_141313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cordon_del_azufre_083D_11651_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cordon_del_azufre_149A_11313_383737 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cordon_del_azufre_149A_11428_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cordon_del_azufre_156D_11344_313232 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cordon_del_azufre_156D_11424_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cordon_del_azufre_156D_11622_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/domuyo_018A_12668_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/domuyo_083D_12636_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/domuyo_091A_12777_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_106A_09048_000909 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_107A_09015_000505 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_107A_09018_000205 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_128D_09016_110500 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/laguna_del_maule_018A_12668_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/laguna_del_maule_083D_12636_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/masaya_055D_07671_161411 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/masaya_157D_07716_131313 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/masaya_165A_07753_101313 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sabancaya_025D_10464_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sabancaya_025D_10613_090705 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sabancaya_047A_10580_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sabancaya_127D_10619_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_miguel_063A_07559_161820 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_miguel_063A_07617_091113 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_miguel_128D_07552_151615 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_miguel_128D_07650_030403 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sangay_018A_09253_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sangay_142D_09147_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sangay_142D_09148_131313 licsalert_sync/south_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_ana_063A_07559_161820 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_ana_063A_07617_091113 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_ana_128D_07531_131313 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_ana_128D_07552_151615 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_ana_136A_07496_151719 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/sierra_negra_106A_09048_000909 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/sierra_negra_106A_09090_000404 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/sierra_negra_128D_09016_110500 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/turrialba_084D_08014_091312 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/turrialba_092A_08084_081213 licsalert_sync/central_america
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_106A_09048_000909 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_107A_09015_000505 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_107A_09018_000205 licsalert_sync/pacific_island
rsync -av --exclude='json' mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_128D_09016_110500 licsalert_sync/pacific_island
