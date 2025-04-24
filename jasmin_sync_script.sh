#!/bin/bash

# Note - this will only work if you can login to Jasmin
# which requires running a command such as:
# eval $(ssh-agent -s); ssh-add ~/.ssh/id_rsa_jasmin
# after you have setup a key pair with Jasmin.


# Trap SIGINT (Ctrl+C)
trap 'echo "Script interrupted, exiting..."; exit 1' INT

echo 'Attempting to rsync alcedo_106A_09048_000909 (0 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/alcedo_106A_09048_000909 licsalert_sync/pacific_island

echo 'Attempting to rsync alcedo_107A_09018_000205 (1 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/alcedo_107A_09018_000205 licsalert_sync/pacific_island

echo 'Attempting to rsync alcedo_128D_09016_110500 (2 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/alcedo_128D_09016_110500 licsalert_sync/pacific_island

echo 'Attempting to rsync askja_009D_02423_131316 (3 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/askja_009D_02423_131316 licsalert_sync/iceland

echo 'Attempting to rsync askja_009D_02504_202119 (4 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/askja_009D_02504_202119 licsalert_sync/iceland

echo 'Attempting to rsync askja_045A_02494_171816 (5 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/askja_045A_02494_171816 licsalert_sync/iceland

echo 'Attempting to rsync askja_045A_02538_121210 (6 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/askja_045A_02538_121210 licsalert_sync/iceland

echo 'Attempting to rsync askja_111D_02489_111313 (7 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/askja_111D_02489_111313 licsalert_sync/iceland

echo 'Attempting to rsync askja_111D_02490_152021 (8 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/askja_111D_02490_152021 licsalert_sync/iceland

echo 'Attempting to rsync askja_118A_02507_211817 (9 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/askja_118A_02507_211817 licsalert_sync/iceland

echo 'Attempting to rsync askja_147A_02466_191712 (10 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/askja_147A_02466_191712 licsalert_sync/iceland

echo 'Attempting to rsync askja_147A_02488_131211 (11 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/askja_147A_02488_131211 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_009D_02423_131316 (12 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_009D_02423_131316 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_009D_02504_202119 (13 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_009D_02504_202119 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_045A_02494_171816 (14 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_045A_02494_171816 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_045A_02538_121210 (15 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_045A_02538_121210 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_111D_02489_111313 (16 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_111D_02489_111313 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_111D_02490_152021 (17 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_111D_02490_152021 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_118A_02507_211817 (18 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_118A_02507_211817 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_147A_02466_191712 (19 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_147A_02466_191712 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_147A_02488_131211 (20 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_147A_02488_131211 licsalert_sync/iceland

echo 'Attempting to rsync campi_flegrei_044A_04913_071213 (21 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/europe/campi_flegrei_044A_04913_071213 licsalert_sync/europe

echo 'Attempting to rsync campi_flegrei_124D_04854_171313 (22 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/europe/campi_flegrei_124D_04854_171313 licsalert_sync/europe

echo 'Attempting to rsync cerro_azul_106A_09048_000909 (23 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/cerro_azul_106A_09048_000909 licsalert_sync/pacific_island

echo 'Attempting to rsync cerro_azul_106A_09090_000404 (24 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/cerro_azul_106A_09090_000404 licsalert_sync/pacific_island

echo 'Attempting to rsync cerro_azul_128D_09016_110500 (25 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/cerro_azul_128D_09016_110500 licsalert_sync/pacific_island

echo 'Attempting to rsync darwin_106A_09048_000909 (26 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/darwin_106A_09048_000909 licsalert_sync/pacific_island

echo 'Attempting to rsync darwin_107A_09015_000505 (27 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/darwin_107A_09015_000505 licsalert_sync/pacific_island

echo 'Attempting to rsync darwin_107A_09018_000205 (28 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/darwin_107A_09018_000205 licsalert_sync/pacific_island

echo 'Attempting to rsync darwin_128D_09016_110500 (29 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/darwin_128D_09016_110500 licsalert_sync/pacific_island

echo 'Attempting to rsync domuyo_018A_12668_131313 (30 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/domuyo_018A_12668_131313 licsalert_sync/south_america

echo 'Attempting to rsync domuyo_083D_12636_131313 (31 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/domuyo_083D_12636_131313 licsalert_sync/south_america

echo 'Attempting to rsync domuyo_091A_12777_131313 (32 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/domuyo_091A_12777_131313 licsalert_sync/south_america

echo 'Attempting to rsync ecuador_106A_09048_000909 (33 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/ecuador_106A_09048_000909 licsalert_sync/pacific_island

echo 'Attempting to rsync ecuador_107A_09015_000505 (34 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/ecuador_107A_09015_000505 licsalert_sync/pacific_island

echo 'Attempting to rsync ecuador_107A_09018_000205 (35 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/ecuador_107A_09018_000205 licsalert_sync/pacific_island

echo 'Attempting to rsync ecuador_128D_09016_110500 (36 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/ecuador_128D_09016_110500 licsalert_sync/pacific_island

echo 'Attempting to rsync etna_044A_05255_141005 (37 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/europe/etna_044A_05255_141005 licsalert_sync/europe

echo 'Attempting to rsync etna_124D_05291_081406 (38 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/europe/etna_124D_05291_081406 licsalert_sync/europe

echo 'Attempting to rsync fernandina_106A_09048_000909 (39 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_106A_09048_000909 licsalert_sync/pacific_island

echo 'Attempting to rsync fernandina_107A_09015_000505 (40 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_107A_09015_000505 licsalert_sync/pacific_island

echo 'Attempting to rsync fernandina_107A_09018_000205 (41 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_107A_09018_000205 licsalert_sync/pacific_island

echo 'Attempting to rsync fernandina_128D_09016_110500 (42 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_128D_09016_110500 licsalert_sync/pacific_island

echo 'Attempting to rsync fujisan_039A_05372_151515 (43 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/fujisan_039A_05372_151515 licsalert_sync/eastern_asia

echo 'Attempting to rsync fujisan_046D_05469_071311 (44 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/fujisan_046D_05469_071311 licsalert_sync/eastern_asia

echo 'Attempting to rsync fujisan_119D_05372_141313 (45 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/fujisan_119D_05372_141313 licsalert_sync/eastern_asia

echo 'Attempting to rsync iztaccihuatl_005A_07021_131313 (46 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/iztaccihuatl_005A_07021_131313 licsalert_sync/north_america

echo 'Attempting to rsync iztaccihuatl_143D_07132_131313 (47 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/iztaccihuatl_143D_07132_131313 licsalert_sync/north_america

echo 'Attempting to rsync iztaccihuatl_143D_07197_212120 (48 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/iztaccihuatl_143D_07197_212120 licsalert_sync/north_america

echo 'Attempting to rsync kilauea_087D_07004_060904 (49 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/kilauea_087D_07004_060904 licsalert_sync/pacific_island

echo 'Attempting to rsync kilauea_124A_06996_091406 (50 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/kilauea_124A_06996_091406 licsalert_sync/pacific_island

echo 'Attempting to rsync krysuvik_016A_02504_162118 (51 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/krysuvik_016A_02504_162118 licsalert_sync/iceland

echo 'Attempting to rsync krysuvik_016A_02562_091313 (52 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/krysuvik_016A_02562_091313 licsalert_sync/iceland

echo 'Attempting to rsync krysuvik_155D_02484_191814 (53 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/krysuvik_155D_02484_191814 licsalert_sync/iceland

echo 'Attempting to rsync krysuvik_155D_02579_100800 (54 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/krysuvik_155D_02579_100800 licsalert_sync/iceland

echo 'Attempting to rsync kverkfjoll_009D_02423_131316 (55 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/kverkfjoll_009D_02423_131316 licsalert_sync/iceland

echo 'Attempting to rsync kverkfjoll_009D_02504_202119 (56 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/kverkfjoll_009D_02504_202119 licsalert_sync/iceland

echo 'Attempting to rsync kverkfjoll_045A_02494_171816 (57 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/kverkfjoll_045A_02494_171816 licsalert_sync/iceland

echo 'Attempting to rsync kverkfjoll_045A_02538_121210 (58 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/kverkfjoll_045A_02538_121210 licsalert_sync/iceland

echo 'Attempting to rsync kverkfjoll_111D_02489_111313 (59 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/kverkfjoll_111D_02489_111313 licsalert_sync/iceland

echo 'Attempting to rsync kverkfjoll_111D_02490_152021 (60 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/kverkfjoll_111D_02490_152021 licsalert_sync/iceland

echo 'Attempting to rsync kverkfjoll_147A_02466_191712 (61 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/kverkfjoll_147A_02466_191712 licsalert_sync/iceland

echo 'Attempting to rsync kverkfjoll_147A_02488_131211 (62 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/kverkfjoll_147A_02488_131211 licsalert_sync/iceland

echo 'Attempting to rsync laguna_del_maule_018A_12668_131313 (63 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/laguna_del_maule_018A_12668_131313 licsalert_sync/south_america

echo 'Attempting to rsync laguna_del_maule_083D_12636_131313 (64 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/laguna_del_maule_083D_12636_131313 licsalert_sync/south_america

echo 'Attempting to rsync mauna_kea_087D_07004_060904 (65 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/mauna_kea_087D_07004_060904 licsalert_sync/pacific_island

echo 'Attempting to rsync mauna_kea_124A_06996_091406 (66 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/mauna_kea_124A_06996_091406 licsalert_sync/pacific_island

echo 'Attempting to rsync mauna_loa_087D_07004_060904 (67 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/mauna_loa_087D_07004_060904 licsalert_sync/pacific_island

echo 'Attempting to rsync mauna_loa_124A_06996_091406 (68 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/mauna_loa_124A_06996_091406 licsalert_sync/pacific_island

echo 'Attempting to rsync pico_002A_05136_020502 (69 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/atlantic_island/pico_002A_05136_020502 licsalert_sync/atlantic_island

echo 'Attempting to rsync pico_082D_05128_030500 (70 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/atlantic_island/pico_082D_05128_030500 licsalert_sync/atlantic_island

echo 'Attempting to rsync reykjanes_016A_02504_162118 (71 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/reykjanes_016A_02504_162118 licsalert_sync/iceland

echo 'Attempting to rsync reykjanes_016A_02562_091313 (72 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/reykjanes_016A_02562_091313 licsalert_sync/iceland

echo 'Attempting to rsync reykjanes_155D_02484_191814 (73 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/reykjanes_155D_02484_191814 licsalert_sync/iceland

echo 'Attempting to rsync reykjanes_155D_02579_100800 (74 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/reykjanes_155D_02579_100800 licsalert_sync/iceland

echo 'Attempting to rsync sierra_negra_106A_09048_000909 (75 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/sierra_negra_106A_09048_000909 licsalert_sync/pacific_island

echo 'Attempting to rsync sierra_negra_106A_09090_000404 (76 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/sierra_negra_106A_09090_000404 licsalert_sync/pacific_island

echo 'Attempting to rsync sierra_negra_128D_09016_110500 (77 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/sierra_negra_128D_09016_110500 licsalert_sync/pacific_island

echo 'Attempting to rsync wolf_106A_09048_000909 (78 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_106A_09048_000909 licsalert_sync/pacific_island

echo 'Attempting to rsync wolf_107A_09015_000505 (79 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_107A_09015_000505 licsalert_sync/pacific_island

echo 'Attempting to rsync wolf_107A_09018_000205 (80 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_107A_09018_000205 licsalert_sync/pacific_island

echo 'Attempting to rsync wolf_128D_09016_110500 (81 of 82)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_128D_09016_110500 licsalert_sync/pacific_island

