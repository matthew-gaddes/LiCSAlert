#!/bin/bash

# Note - this will only work if you can login to Jasmin
# which requires running a command such as:
# eval $(ssh-agent -s); ssh-add ~/.ssh/id_rsa_jasmin
# after you have setup a key pair with Jasmin.


# Trap SIGINT (Ctrl+C)
trap 'echo "Script interrupted, exiting..."; exit 1' INT

echo 'Attempting to rsync stromboli_044A_05129_010307 (0 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/europe/stromboli_044A_05129_010307 licsalert_sync/europe

echo 'Attempting to rsync stromboli_124D_05178_010601 (1 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/europe/stromboli_124D_05178_010601 licsalert_sync/europe

echo 'Attempting to rsync etna_044A_05255_141005 (2 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/europe/etna_044A_05255_141005 licsalert_sync/europe

echo 'Attempting to rsync etna_124D_05291_081406 (3 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/europe/etna_124D_05291_081406 licsalert_sync/europe

echo 'Attempting to rsync erta_ale_014A_07688_131313 (4 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/erta_ale_014A_07688_131313 licsalert_sync/africa

echo 'Attempting to rsync erta_ale_079D_07694_131313 (5 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/erta_ale_079D_07694_131313 licsalert_sync/africa

echo 'Attempting to rsync dabbahu_014A_07688_131313 (6 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/dabbahu_014A_07688_131313 licsalert_sync/africa

echo 'Attempting to rsync dabbahu_079D_07694_131313 (7 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/dabbahu_079D_07694_131313 licsalert_sync/africa

echo 'Attempting to rsync manda_hararo_014A_07688_131313 (8 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/manda_hararo_014A_07688_131313 licsalert_sync/africa

echo 'Attempting to rsync manda_hararo_014A_07884_121313 (9 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/manda_hararo_014A_07884_121313 licsalert_sync/africa

echo 'Attempting to rsync manda_hararo_079D_07694_131313 (10 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/manda_hararo_079D_07694_131313 licsalert_sync/africa

echo 'Attempting to rsync ol_doinyo_lengai_130A_09212_131313 (11 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/ol_doinyo_lengai_130A_09212_131313 licsalert_sync/africa

echo 'Attempting to rsync ol_doinyo_lengai_152D_09315_131313 (12 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/ol_doinyo_lengai_152D_09315_131313 licsalert_sync/africa

echo 'Attempting to rsync nyamuragira_021D_09150_131313 (13 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/nyamuragira_021D_09150_131313 licsalert_sync/africa

echo 'Attempting to rsync nyamuragira_174A_09133_131313 (14 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/nyamuragira_174A_09133_131313 licsalert_sync/africa

echo 'Attempting to rsync nyiragongo_021D_09150_131313 (15 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/nyiragongo_021D_09150_131313 licsalert_sync/africa

echo 'Attempting to rsync nyiragongo_174A_09133_131313 (16 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/africa/nyiragongo_174A_09133_131313 licsalert_sync/africa

echo 'Attempting to rsync piton_de_la_fournaise_144A_11114_000304 (17 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/indian_island/piton_de_la_fournaise_144A_11114_000304 licsalert_sync/indian_island

echo 'Attempting to rsync piton_de_la_fournaise_151D_11112_000300 (18 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/indian_island/piton_de_la_fournaise_151D_11112_000300 licsalert_sync/indian_island

echo 'Attempting to rsync white_island_008A_12761_090302 (19 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/oceania/white_island_008A_12761_090302 licsalert_sync/oceania

echo 'Attempting to rsync langila_162D_09619_091314 (20 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/oceania/langila_162D_09619_091314 licsalert_sync/oceania

echo 'Attempting to rsync bagana_096A_09625_110000 (21 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/oceania/bagana_096A_09625_110000 licsalert_sync/oceania

echo 'Attempting to rsync sinabung_062D_08629_071113 (22 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/sinabung_062D_08629_071113 licsalert_sync/southeast_asia

echo 'Attempting to rsync sinabung_143A_08674_131308 (23 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/sinabung_143A_08674_131308 licsalert_sync/southeast_asia

echo 'Attempting to rsync marapi_091D_09086_131313 (24 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/marapi_091D_09086_131313 licsalert_sync/southeast_asia

echo 'Attempting to rsync marapi_164D_08941_131113 (25 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/marapi_164D_08941_131113 licsalert_sync/southeast_asia

echo 'Attempting to rsync talang_091D_09086_131313 (26 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/talang_091D_09086_131313 licsalert_sync/southeast_asia

echo 'Attempting to rsync kaba_018D_09287_131313 (27 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/kaba_018D_09287_131313 licsalert_sync/southeast_asia

echo 'Attempting to rsync kaba_069A_09457_021206 (28 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/kaba_069A_09457_021206 licsalert_sync/southeast_asia

echo 'Attempting to rsync dempo_018D_09415_060403 (29 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/dempo_018D_09415_060403 licsalert_sync/southeast_asia

echo 'Attempting to rsync dempo_069A_09283_080913 (30 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/dempo_069A_09283_080913 licsalert_sync/southeast_asia

echo 'Attempting to rsync dempo_069A_09457_021206 (31 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/dempo_069A_09457_021206 licsalert_sync/southeast_asia

echo 'Attempting to rsync tangkubanparahu_098A_09673_121312 (32 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/tangkubanparahu_098A_09673_121312 licsalert_sync/southeast_asia

echo 'Attempting to rsync tangkubanparahu_149D_09700_081210 (33 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/tangkubanparahu_149D_09700_081210 licsalert_sync/southeast_asia

echo 'Attempting to rsync papandayan_098A_09673_121312 (34 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/papandayan_098A_09673_121312 licsalert_sync/southeast_asia

echo 'Attempting to rsync papandayan_149D_09700_081210 (35 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/papandayan_149D_09700_081210 licsalert_sync/southeast_asia

echo 'Attempting to rsync slamet_025A_09718_120810 (36 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/slamet_025A_09718_120810 licsalert_sync/southeast_asia

echo 'Attempting to rsync slamet_076D_09725_121107 (37 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/slamet_076D_09725_121107 licsalert_sync/southeast_asia

echo 'Attempting to rsync slamet_149D_09700_081210 (38 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/slamet_149D_09700_081210 licsalert_sync/southeast_asia

echo 'Attempting to rsync dieng_volcanic_complex_025A_09718_120810 (39 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/dieng_volcanic_complex_025A_09718_120810 licsalert_sync/southeast_asia

echo 'Attempting to rsync dieng_volcanic_complex_076D_09725_121107 (40 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/dieng_volcanic_complex_076D_09725_121107 licsalert_sync/southeast_asia

echo 'Attempting to rsync merapi_076D_09725_121107 (41 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/merapi_076D_09725_121107 licsalert_sync/southeast_asia

echo 'Attempting to rsync merapi_127A_09749_121312 (42 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/merapi_127A_09749_121312 licsalert_sync/southeast_asia

echo 'Attempting to rsync kelut_003D_09757_111111 (43 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/kelut_003D_09757_111111 licsalert_sync/southeast_asia

echo 'Attempting to rsync kelut_054A_09773_111213 (44 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/kelut_054A_09773_111213 licsalert_sync/southeast_asia

echo 'Attempting to rsync kelut_127A_09749_121312 (45 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/kelut_127A_09749_121312 licsalert_sync/southeast_asia

echo 'Attempting to rsync semeru_003D_09757_111111 (46 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/semeru_003D_09757_111111 licsalert_sync/southeast_asia

echo 'Attempting to rsync semeru_054A_09773_111213 (47 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/semeru_054A_09773_111213 licsalert_sync/southeast_asia

echo 'Attempting to rsync tengger_caldera_003D_09757_111111 (48 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/tengger_caldera_003D_09757_111111 licsalert_sync/southeast_asia

echo 'Attempting to rsync tengger_caldera_054A_09773_111213 (49 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/tengger_caldera_054A_09773_111213 licsalert_sync/southeast_asia

echo 'Attempting to rsync raung_054A_09773_111213 (50 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/raung_054A_09773_111213 licsalert_sync/southeast_asia

echo 'Attempting to rsync raung_105D_09782_131111 (51 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/raung_105D_09782_131111 licsalert_sync/southeast_asia

echo 'Attempting to rsync ijen_054A_09773_111213 (52 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/ijen_054A_09773_111213 licsalert_sync/southeast_asia

echo 'Attempting to rsync ijen_105D_09782_131111 (53 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/ijen_105D_09782_131111 licsalert_sync/southeast_asia

echo 'Attempting to rsync batur_032D_09854_070505 (54 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/batur_032D_09854_070505 licsalert_sync/southeast_asia

echo 'Attempting to rsync batur_105D_09782_131111 (55 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/batur_105D_09782_131111 licsalert_sync/southeast_asia

echo 'Attempting to rsync batur_156A_09814_081406 (56 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/batur_156A_09814_081406 licsalert_sync/southeast_asia

echo 'Attempting to rsync agung_032D_09854_070505 (57 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/agung_032D_09854_070505 licsalert_sync/southeast_asia

echo 'Attempting to rsync agung_156A_09814_081406 (58 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/agung_156A_09814_081406 licsalert_sync/southeast_asia

echo 'Attempting to rsync rinjani_032D_09854_070505 (59 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/rinjani_032D_09854_070505 licsalert_sync/southeast_asia

echo 'Attempting to rsync rinjani_156A_09814_081406 (60 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/rinjani_156A_09814_081406 licsalert_sync/southeast_asia

echo 'Attempting to rsync sangeang_api_010A_09915_111413 (61 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/sangeang_api_010A_09915_111413 licsalert_sync/southeast_asia

echo 'Attempting to rsync sangeang_api_061D_09920_121210 (62 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/sangeang_api_061D_09920_121210 licsalert_sync/southeast_asia

echo 'Attempting to rsync ranakah_010A_09915_111413 (63 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/ranakah_010A_09915_111413 licsalert_sync/southeast_asia

echo 'Attempting to rsync ranakah_061D_09920_121210 (64 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/ranakah_061D_09920_121210 licsalert_sync/southeast_asia

echo 'Attempting to rsync inielika_010A_09915_111413 (65 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/inielika_010A_09915_111413 licsalert_sync/southeast_asia

echo 'Attempting to rsync inielika_061D_09920_121210 (66 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/inielika_061D_09920_121210 licsalert_sync/southeast_asia

echo 'Attempting to rsync inielika_112A_09831_050508 (67 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/inielika_112A_09831_050508 licsalert_sync/southeast_asia

echo 'Attempting to rsync paluweh_112A_09831_050508 (68 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/paluweh_112A_09831_050508 licsalert_sync/southeast_asia

echo 'Attempting to rsync paluweh_163D_09852_050706 (69 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/paluweh_163D_09852_050706 licsalert_sync/southeast_asia

echo 'Attempting to rsync egon_112A_09831_050508 (70 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/egon_112A_09831_050508 licsalert_sync/southeast_asia

echo 'Attempting to rsync egon_163D_09852_050706 (71 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/egon_163D_09852_050706 licsalert_sync/southeast_asia

echo 'Attempting to rsync lewotobi_163D_09852_050706 (72 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/lewotobi_163D_09852_050706 licsalert_sync/southeast_asia

echo 'Attempting to rsync leroboleng_039A_09827_040302 (73 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/leroboleng_039A_09827_040302 licsalert_sync/southeast_asia

echo 'Attempting to rsync leroboleng_163D_09852_050706 (74 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/leroboleng_163D_09852_050706 licsalert_sync/southeast_asia

echo 'Attempting to rsync iliboleng_039A_09827_040302 (75 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/iliboleng_039A_09827_040302 licsalert_sync/southeast_asia

echo 'Attempting to rsync iliboleng_039A_09931_061510 (76 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/iliboleng_039A_09931_061510 licsalert_sync/southeast_asia

echo 'Attempting to rsync iliboleng_090D_09898_131213 (77 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/iliboleng_090D_09898_131213 licsalert_sync/southeast_asia

echo 'Attempting to rsync iliboleng_163D_09852_050706 (78 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/iliboleng_163D_09852_050706 licsalert_sync/southeast_asia

echo 'Attempting to rsync lewotolo_039A_09827_040302 (79 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/lewotolo_039A_09827_040302 licsalert_sync/southeast_asia

echo 'Attempting to rsync lewotolo_039A_09931_061510 (80 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/lewotolo_039A_09931_061510 licsalert_sync/southeast_asia

echo 'Attempting to rsync lewotolo_090D_09898_131213 (81 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/lewotolo_090D_09898_131213 licsalert_sync/southeast_asia

echo 'Attempting to rsync iliwerung_039A_09827_040302 (82 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/iliwerung_039A_09827_040302 licsalert_sync/southeast_asia

echo 'Attempting to rsync iliwerung_039A_09931_061510 (83 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/iliwerung_039A_09931_061510 licsalert_sync/southeast_asia

echo 'Attempting to rsync iliwerung_090D_09898_131213 (84 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/iliwerung_090D_09898_131213 licsalert_sync/southeast_asia

echo 'Attempting to rsync sirung_039A_09827_040302 (85 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/sirung_039A_09827_040302 licsalert_sync/southeast_asia

echo 'Attempting to rsync sirung_039A_09931_061510 (86 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/sirung_039A_09931_061510 licsalert_sync/southeast_asia

echo 'Attempting to rsync sirung_090D_09898_131213 (87 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/sirung_090D_09898_131213 licsalert_sync/southeast_asia

echo 'Attempting to rsync soputan_163D_08898_091005 (88 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/soputan_163D_08898_091005 licsalert_sync/southeast_asia

echo 'Attempting to rsync lokon-empung_163D_08898_091005 (89 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/lokon-empung_163D_08898_091005 licsalert_sync/southeast_asia

echo 'Attempting to rsync kanlaon_069A_07955_090909 (90 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/kanlaon_069A_07955_090909 licsalert_sync/southeast_asia

echo 'Attempting to rsync kanlaon_134D_07892_131210 (91 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/kanlaon_134D_07892_131210 licsalert_sync/southeast_asia

echo 'Attempting to rsync bulusan_069A_07674_131310 (92 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/bulusan_069A_07674_131310 licsalert_sync/southeast_asia

echo 'Attempting to rsync mayon_069A_07674_131310 (93 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/mayon_069A_07674_131310 licsalert_sync/southeast_asia

echo 'Attempting to rsync mayon_134D_07688_131313 (94 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/mayon_134D_07688_131313 licsalert_sync/southeast_asia

echo 'Attempting to rsync pinatubo_032D_07536_111313 (95 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/pinatubo_032D_07536_111313 licsalert_sync/southeast_asia

echo 'Attempting to rsync pinatubo_142A_07411_111313 (96 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/southeast_asia/pinatubo_142A_07411_111313 licsalert_sync/southeast_asia

echo 'Attempting to rsync suwanosejima_163D_05934_091414 (97 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/suwanosejima_163D_05934_091414 licsalert_sync/eastern_asia

echo 'Attempting to rsync kuchinoerabujima_054A_05884_001412 (98 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/kuchinoerabujima_054A_05884_001412 licsalert_sync/eastern_asia

echo 'Attempting to rsync kuchinoerabujima_163D_05934_091414 (99 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/kuchinoerabujima_163D_05934_091414 licsalert_sync/eastern_asia

echo 'Attempting to rsync aira_054A_05884_001412 (100 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/aira_054A_05884_001412 licsalert_sync/eastern_asia

echo 'Attempting to rsync aira_156A_05796_080500 (101 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/aira_156A_05796_080500 licsalert_sync/eastern_asia

echo 'Attempting to rsync aira_163D_05736_131313 (102 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/aira_163D_05736_131313 licsalert_sync/eastern_asia

echo 'Attempting to rsync aira_163D_05934_091414 (103 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/aira_163D_05934_091414 licsalert_sync/eastern_asia

echo 'Attempting to rsync kirishimayama_054A_05884_001412 (104 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/kirishimayama_054A_05884_001412 licsalert_sync/eastern_asia

echo 'Attempting to rsync kirishimayama_156A_05796_080500 (105 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/kirishimayama_156A_05796_080500 licsalert_sync/eastern_asia

echo 'Attempting to rsync kirishimayama_163D_05736_131313 (106 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/kirishimayama_163D_05736_131313 licsalert_sync/eastern_asia

echo 'Attempting to rsync unzendake_054A_05688_051313 (107 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/unzendake_054A_05688_051313 licsalert_sync/eastern_asia

echo 'Attempting to rsync unzendake_156A_05639_121313 (108 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/unzendake_156A_05639_121313 licsalert_sync/eastern_asia

echo 'Attempting to rsync unzendake_163D_05736_131313 (109 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/unzendake_163D_05736_131313 licsalert_sync/eastern_asia

echo 'Attempting to rsync asosan_156A_05639_121313 (110 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/asosan_156A_05639_121313 licsalert_sync/eastern_asia

echo 'Attempting to rsync asosan_163D_05736_131313 (111 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/asosan_163D_05736_131313 licsalert_sync/eastern_asia

echo 'Attempting to rsync kujusan_156A_05639_121313 (112 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/kujusan_156A_05639_121313 licsalert_sync/eastern_asia

echo 'Attempting to rsync kujusan_163D_05736_131313 (113 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/kujusan_163D_05736_131313 licsalert_sync/eastern_asia

echo 'Attempting to rsync hakoneyama_039A_05372_151515 (114 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/hakoneyama_039A_05372_151515 licsalert_sync/eastern_asia

echo 'Attempting to rsync hakoneyama_046D_05469_071311 (115 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/hakoneyama_046D_05469_071311 licsalert_sync/eastern_asia

echo 'Attempting to rsync hakoneyama_119D_05372_141313 (116 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/hakoneyama_119D_05372_141313 licsalert_sync/eastern_asia

echo 'Attempting to rsync ontakesan_112A_05454_131213 (117 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/ontakesan_112A_05454_131213 licsalert_sync/eastern_asia

echo 'Attempting to rsync ontakesan_119D_05372_141313 (118 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/ontakesan_119D_05372_141313 licsalert_sync/eastern_asia

echo 'Attempting to rsync yakedake_112A_05454_131213 (119 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/yakedake_112A_05454_131213 licsalert_sync/eastern_asia

echo 'Attempting to rsync yakedake_119D_05372_141313 (120 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/yakedake_119D_05372_141313 licsalert_sync/eastern_asia

echo 'Attempting to rsync niigata-yakeyama_039A_05372_151515 (121 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/niigata-yakeyama_039A_05372_151515 licsalert_sync/eastern_asia

echo 'Attempting to rsync niigata-yakeyama_112A_05300_030808 (122 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/niigata-yakeyama_112A_05300_030808 licsalert_sync/eastern_asia

echo 'Attempting to rsync niigata-yakeyama_119D_05372_141313 (123 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/niigata-yakeyama_119D_05372_141313 licsalert_sync/eastern_asia

echo 'Attempting to rsync asamayama_039A_05372_151515 (124 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/asamayama_039A_05372_151515 licsalert_sync/eastern_asia

echo 'Attempting to rsync asamayama_112A_05454_131213 (125 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/asamayama_112A_05454_131213 licsalert_sync/eastern_asia

echo 'Attempting to rsync asamayama_119D_05372_141313 (126 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/asamayama_119D_05372_141313 licsalert_sync/eastern_asia

echo 'Attempting to rsync kusatsu-shiranesan_039A_05372_151515 (127 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/kusatsu-shiranesan_039A_05372_151515 licsalert_sync/eastern_asia

echo 'Attempting to rsync kusatsu-shiranesan_119D_05372_141313 (128 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/kusatsu-shiranesan_119D_05372_141313 licsalert_sync/eastern_asia

echo 'Attempting to rsync adatarayama_046D_05292_131313 (129 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/adatarayama_046D_05292_131313 licsalert_sync/eastern_asia

echo 'Attempting to rsync adatarayama_141A_05337_160700 (130 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/adatarayama_141A_05337_160700 licsalert_sync/eastern_asia

echo 'Attempting to rsync akita-yakeyama_046D_04907_081313 (131 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/akita-yakeyama_046D_04907_081313 licsalert_sync/eastern_asia

echo 'Attempting to rsync akita-yakeyama_046D_05096_121313 (132 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/akita-yakeyama_046D_05096_121313 licsalert_sync/eastern_asia

echo 'Attempting to rsync akita-yakeyama_141A_04925_131313 (133 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/akita-yakeyama_141A_04925_131313 licsalert_sync/eastern_asia

echo 'Attempting to rsync izu-oshima_039A_05599_101702 (134 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/izu-oshima_039A_05599_101702 licsalert_sync/eastern_asia

echo 'Attempting to rsync izu-oshima_046D_05469_071311 (135 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/izu-oshima_046D_05469_071311 licsalert_sync/eastern_asia

echo 'Attempting to rsync miyakejima_039A_05599_101702 (136 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/miyakejima_039A_05599_101702 licsalert_sync/eastern_asia

echo 'Attempting to rsync miyakejima_046D_05469_071311 (137 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/miyakejima_046D_05469_071311 licsalert_sync/eastern_asia

echo 'Attempting to rsync hokkaido-komagatake_046D_04690_121308 (138 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/hokkaido-komagatake_046D_04690_121308 licsalert_sync/eastern_asia

echo 'Attempting to rsync hokkaido-komagatake_046D_04907_081313 (139 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/hokkaido-komagatake_046D_04907_081313 licsalert_sync/eastern_asia

echo 'Attempting to rsync hokkaido-komagatake_068A_04799_130306 (140 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/hokkaido-komagatake_068A_04799_130306 licsalert_sync/eastern_asia

echo 'Attempting to rsync hokkaido-komagatake_141A_04773_020709 (141 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/hokkaido-komagatake_141A_04773_020709 licsalert_sync/eastern_asia

echo 'Attempting to rsync toya_046D_04690_121308 (142 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/toya_046D_04690_121308 licsalert_sync/eastern_asia

echo 'Attempting to rsync toya_068A_04667_071312 (143 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/toya_068A_04667_071312 licsalert_sync/eastern_asia

echo 'Attempting to rsync toya_068A_04799_130306 (144 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/toya_068A_04799_130306 licsalert_sync/eastern_asia

echo 'Attempting to rsync toya_141A_04773_020709 (145 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/toya_141A_04773_020709 licsalert_sync/eastern_asia

echo 'Attempting to rsync tokachidake_046D_04690_121308 (146 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/tokachidake_046D_04690_121308 licsalert_sync/eastern_asia

echo 'Attempting to rsync tokachidake_068A_04667_071312 (147 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/tokachidake_068A_04667_071312 licsalert_sync/eastern_asia

echo 'Attempting to rsync tokachidake_148D_04631_101009 (148 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/tokachidake_148D_04631_101009 licsalert_sync/eastern_asia

echo 'Attempting to rsync tokachidake_170A_04675_131008 (149 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/tokachidake_170A_04675_131008 licsalert_sync/eastern_asia

echo 'Attempting to rsync akan_148D_04618_101013 (150 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/akan_148D_04618_101013 licsalert_sync/eastern_asia

echo 'Attempting to rsync akan_148D_04631_101009 (151 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/akan_148D_04631_101009 licsalert_sync/eastern_asia

echo 'Attempting to rsync akan_170A_04675_131008 (152 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/eastern_asia/akan_170A_04675_131008 licsalert_sync/eastern_asia

echo 'Attempting to rsync chirpoi_053A_04342_030304 (153 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/chirpoi_053A_04342_030304 licsalert_sync/northern_asia

echo 'Attempting to rsync kambalny_111A_03775_131006 (154 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/kambalny_111A_03775_131006 licsalert_sync/northern_asia

echo 'Attempting to rsync koryaksky_038A_03576_131307 (155 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/koryaksky_038A_03576_131307 licsalert_sync/northern_asia

echo 'Attempting to rsync koryaksky_060D_03564_141313 (156 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/koryaksky_060D_03564_141313 licsalert_sync/northern_asia

echo 'Attempting to rsync koryaksky_060D_03718_021307 (157 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/koryaksky_060D_03718_021307 licsalert_sync/northern_asia

echo 'Attempting to rsync koryaksky_111A_03605_131313 (158 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/koryaksky_111A_03605_131313 licsalert_sync/northern_asia

echo 'Attempting to rsync avachinsky_038A_03576_131307 (159 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/avachinsky_038A_03576_131307 licsalert_sync/northern_asia

echo 'Attempting to rsync avachinsky_060D_03564_141313 (160 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/avachinsky_060D_03564_141313 licsalert_sync/northern_asia

echo 'Attempting to rsync avachinsky_060D_03718_021307 (161 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/avachinsky_060D_03718_021307 licsalert_sync/northern_asia

echo 'Attempting to rsync avachinsky_111A_03605_131313 (162 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/avachinsky_111A_03605_131313 licsalert_sync/northern_asia

echo 'Attempting to rsync zhupanovsky_038A_03576_131307 (163 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/zhupanovsky_038A_03576_131307 licsalert_sync/northern_asia

echo 'Attempting to rsync zhupanovsky_060D_03564_141313 (164 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/zhupanovsky_060D_03564_141313 licsalert_sync/northern_asia

echo 'Attempting to rsync zhupanovsky_111A_03605_131313 (165 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/zhupanovsky_111A_03605_131313 licsalert_sync/northern_asia

echo 'Attempting to rsync karymsky_038A_03576_131307 (166 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/karymsky_038A_03576_131307 licsalert_sync/northern_asia

echo 'Attempting to rsync karymsky_060D_03564_141313 (167 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/karymsky_060D_03564_141313 licsalert_sync/northern_asia

echo 'Attempting to rsync bezymianny_038A_03435_080808 (168 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/bezymianny_038A_03435_080808 licsalert_sync/northern_asia

echo 'Attempting to rsync bezymianny_060D_03363_131313 (169 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/bezymianny_060D_03363_131313 licsalert_sync/northern_asia

echo 'Attempting to rsync bezymianny_140A_03370_131312 (170 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/bezymianny_140A_03370_131312 licsalert_sync/northern_asia

echo 'Attempting to rsync klyuchevskoy_038A_03435_080808 (171 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/klyuchevskoy_038A_03435_080808 licsalert_sync/northern_asia

echo 'Attempting to rsync klyuchevskoy_060D_03363_131313 (172 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/klyuchevskoy_060D_03363_131313 licsalert_sync/northern_asia

echo 'Attempting to rsync klyuchevskoy_140A_03370_131312 (173 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/klyuchevskoy_140A_03370_131312 licsalert_sync/northern_asia

echo 'Attempting to rsync sheveluch_060D_03363_131313 (174 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/sheveluch_060D_03363_131313 licsalert_sync/northern_asia

echo 'Attempting to rsync sheveluch_140A_03370_131312 (175 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/northern_asia/sheveluch_140A_03370_131312 licsalert_sync/northern_asia

echo 'Attempting to rsync great_sitkin_037A_03794_030304 (176 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/great_sitkin_037A_03794_030304 licsalert_sync/north_america

echo 'Attempting to rsync great_sitkin_059D_03819_050404 (177 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/great_sitkin_059D_03819_050404 licsalert_sync/north_america

echo 'Attempting to rsync great_sitkin_110A_03818_040303 (178 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/great_sitkin_110A_03818_040303 licsalert_sync/north_america

echo 'Attempting to rsync great_sitkin_161D_03795_040404 (179 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/great_sitkin_161D_03795_040404 licsalert_sync/north_america

echo 'Attempting to rsync cleveland_066A_03735_010103 (180 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/cleveland_066A_03735_010103 licsalert_sync/north_america

echo 'Attempting to rsync cleveland_066A_03750_030203 (181 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/cleveland_066A_03750_030203 licsalert_sync/north_america

echo 'Attempting to rsync cleveland_117D_03705_010404 (182 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/cleveland_117D_03705_010404 licsalert_sync/north_america

echo 'Attempting to rsync bogoslof_044D_03600_070805 (183 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/bogoslof_044D_03600_070805 licsalert_sync/north_america

echo 'Attempting to rsync bogoslof_044D_03605_070807 (184 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/bogoslof_044D_03605_070807 licsalert_sync/north_america

echo 'Attempting to rsync bogoslof_095A_03649_040505 (185 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/bogoslof_095A_03649_040505 licsalert_sync/north_america

echo 'Attempting to rsync shishaldin_124A_03521_060410 (186 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/shishaldin_124A_03521_060410 licsalert_sync/north_america

echo 'Attempting to rsync pavlof_073D_03460_121207 (187 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/pavlof_073D_03460_121207 licsalert_sync/north_america

echo 'Attempting to rsync pavlof_073D_03461_000300 (188 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/pavlof_073D_03461_000300 licsalert_sync/north_america

echo 'Attempting to rsync pavlof_124A_03521_060410 (189 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/pavlof_124A_03521_060410 licsalert_sync/north_america

echo 'Attempting to rsync pavlof_153A_03412_061209 (190 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/pavlof_153A_03412_061209 licsalert_sync/north_america

echo 'Attempting to rsync veniaminof_102D_03341_121118 (191 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/veniaminof_102D_03341_121118 licsalert_sync/north_america

echo 'Attempting to rsync veniaminof_153A_03412_061209 (192 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/veniaminof_153A_03412_061209 licsalert_sync/north_america

echo 'Attempting to rsync st_helens_137A_04439_131313 (193 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/st_helens_137A_04439_131313 licsalert_sync/north_america

echo 'Attempting to rsync kilauea_087D_07004_060904 (194 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/kilauea_087D_07004_060904 licsalert_sync/pacific_island

echo 'Attempting to rsync kilauea_124A_06996_091406 (195 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/kilauea_124A_06996_091406 licsalert_sync/pacific_island

echo 'Attempting to rsync colima_012D_06936_131313 (196 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/colima_012D_06936_131313 licsalert_sync/north_america

echo 'Attempting to rsync colima_012D_06978_201716 (197 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/colima_012D_06978_201716 licsalert_sync/north_america

echo 'Attempting to rsync colima_012D_07075_080504 (198 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/colima_012D_07075_080504 licsalert_sync/north_america

echo 'Attempting to rsync colima_049A_07091_051113 (199 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/colima_049A_07091_051113 licsalert_sync/north_america

echo 'Attempting to rsync popocatepetl_005A_07021_131313 (200 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/popocatepetl_005A_07021_131313 licsalert_sync/north_america

echo 'Attempting to rsync popocatepetl_143D_07132_131313 (201 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/popocatepetl_143D_07132_131313 licsalert_sync/north_america

echo 'Attempting to rsync popocatepetl_143D_07197_212120 (202 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/north_america/popocatepetl_143D_07197_212120 licsalert_sync/north_america

echo 'Attempting to rsync santa_maria_026D_07514_131212 (203 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_maria_026D_07514_131212 licsalert_sync/central_america

echo 'Attempting to rsync santa_maria_026D_07515_131212 (204 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_maria_026D_07515_131212 licsalert_sync/central_america

echo 'Attempting to rsync santa_maria_099D_07414_161310 (205 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_maria_099D_07414_161310 licsalert_sync/central_america

echo 'Attempting to rsync santa_maria_136A_07496_151719 (206 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_maria_136A_07496_151719 licsalert_sync/central_america

echo 'Attempting to rsync santa_maria_136A_07497_151719 (207 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_maria_136A_07497_151719 licsalert_sync/central_america

echo 'Attempting to rsync santa_maria_136A_07522_091613 (208 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_maria_136A_07522_091613 licsalert_sync/central_america

echo 'Attempting to rsync fuego_026D_07514_131212 (209 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/fuego_026D_07514_131212 licsalert_sync/central_america

echo 'Attempting to rsync fuego_026D_07515_131212 (210 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/fuego_026D_07515_131212 licsalert_sync/central_america

echo 'Attempting to rsync fuego_136A_07497_151719 (211 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/fuego_136A_07497_151719 licsalert_sync/central_america

echo 'Attempting to rsync fuego_136A_07522_091613 (212 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/fuego_136A_07522_091613 licsalert_sync/central_america

echo 'Attempting to rsync pacaya_026D_07514_131212 (213 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/pacaya_026D_07514_131212 licsalert_sync/central_america

echo 'Attempting to rsync pacaya_026D_07515_131212 (214 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/pacaya_026D_07515_131212 licsalert_sync/central_america

echo 'Attempting to rsync pacaya_136A_07496_151719 (215 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/pacaya_136A_07496_151719 licsalert_sync/central_america

echo 'Attempting to rsync pacaya_136A_07497_151719 (216 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/pacaya_136A_07497_151719 licsalert_sync/central_america

echo 'Attempting to rsync pacaya_136A_07522_091613 (217 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/pacaya_136A_07522_091613 licsalert_sync/central_america

echo 'Attempting to rsync santa_ana_063A_07559_161820 (218 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_ana_063A_07559_161820 licsalert_sync/central_america

echo 'Attempting to rsync santa_ana_063A_07617_091113 (219 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_ana_063A_07617_091113 licsalert_sync/central_america

echo 'Attempting to rsync santa_ana_128D_07531_131313 (220 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_ana_128D_07531_131313 licsalert_sync/central_america

echo 'Attempting to rsync santa_ana_128D_07552_151615 (221 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_ana_128D_07552_151615 licsalert_sync/central_america

echo 'Attempting to rsync santa_ana_136A_07496_151719 (222 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/santa_ana_136A_07496_151719 licsalert_sync/central_america

echo 'Attempting to rsync san_salvador_063A_07559_161820 (223 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_salvador_063A_07559_161820 licsalert_sync/central_america

echo 'Attempting to rsync san_salvador_063A_07617_091113 (224 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_salvador_063A_07617_091113 licsalert_sync/central_america

echo 'Attempting to rsync san_salvador_128D_07531_131313 (225 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_salvador_128D_07531_131313 licsalert_sync/central_america

echo 'Attempting to rsync san_salvador_128D_07552_151615 (226 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_salvador_128D_07552_151615 licsalert_sync/central_america

echo 'Attempting to rsync san_salvador_128D_07650_030403 (227 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_salvador_128D_07650_030403 licsalert_sync/central_america

echo 'Attempting to rsync san_miguel_063A_07559_161820 (228 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_miguel_063A_07559_161820 licsalert_sync/central_america

echo 'Attempting to rsync san_miguel_063A_07617_091113 (229 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_miguel_063A_07617_091113 licsalert_sync/central_america

echo 'Attempting to rsync san_miguel_128D_07552_151615 (230 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_miguel_128D_07552_151615 licsalert_sync/central_america

echo 'Attempting to rsync san_miguel_128D_07650_030403 (231 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_miguel_128D_07650_030403 licsalert_sync/central_america

echo 'Attempting to rsync san_cristobal_055D_07654_131311 (232 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_cristobal_055D_07654_131311 licsalert_sync/central_america

echo 'Attempting to rsync san_cristobal_055D_07671_161411 (233 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_cristobal_055D_07671_161411 licsalert_sync/central_america

echo 'Attempting to rsync san_cristobal_165A_07753_101313 (234 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/san_cristobal_165A_07753_101313 licsalert_sync/central_america

echo 'Attempting to rsync telica_055D_07654_131311 (235 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/telica_055D_07654_131311 licsalert_sync/central_america

echo 'Attempting to rsync telica_055D_07671_161411 (236 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/telica_055D_07671_161411 licsalert_sync/central_america

echo 'Attempting to rsync telica_165A_07753_101313 (237 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/telica_165A_07753_101313 licsalert_sync/central_america

echo 'Attempting to rsync cerro_negro_055D_07654_131311 (238 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/cerro_negro_055D_07654_131311 licsalert_sync/central_america

echo 'Attempting to rsync cerro_negro_055D_07671_161411 (239 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/cerro_negro_055D_07671_161411 licsalert_sync/central_america

echo 'Attempting to rsync cerro_negro_165A_07753_101313 (240 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/cerro_negro_165A_07753_101313 licsalert_sync/central_america

echo 'Attempting to rsync momotombo_055D_07654_131311 (241 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/momotombo_055D_07654_131311 licsalert_sync/central_america

echo 'Attempting to rsync momotombo_055D_07671_161411 (242 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/momotombo_055D_07671_161411 licsalert_sync/central_america

echo 'Attempting to rsync momotombo_165A_07753_101313 (243 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/momotombo_165A_07753_101313 licsalert_sync/central_america

echo 'Attempting to rsync masaya_055D_07671_161411 (244 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/masaya_055D_07671_161411 licsalert_sync/central_america

echo 'Attempting to rsync masaya_157D_07716_131313 (245 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/masaya_157D_07716_131313 licsalert_sync/central_america

echo 'Attempting to rsync masaya_165A_07753_101313 (246 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/masaya_165A_07753_101313 licsalert_sync/central_america

echo 'Attempting to rsync concepcion_157D_07909_131307 (247 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/concepcion_157D_07909_131307 licsalert_sync/central_america

echo 'Attempting to rsync rincon_de_la_vieja_157D_07909_131307 (248 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/rincon_de_la_vieja_157D_07909_131307 licsalert_sync/central_america

echo 'Attempting to rsync rincon_de_la_vieja_165A_07945_001113 (249 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/rincon_de_la_vieja_165A_07945_001113 licsalert_sync/central_america

echo 'Attempting to rsync poas_084D_08014_091312 (250 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/poas_084D_08014_091312 licsalert_sync/central_america

echo 'Attempting to rsync poas_084D_08033_121612 (251 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/poas_084D_08033_121612 licsalert_sync/central_america

echo 'Attempting to rsync poas_092A_08036_141719 (252 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/poas_092A_08036_141719 licsalert_sync/central_america

echo 'Attempting to rsync irazu_084D_08014_091312 (253 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/irazu_084D_08014_091312 licsalert_sync/central_america

echo 'Attempting to rsync irazu_084D_08033_121612 (254 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/irazu_084D_08033_121612 licsalert_sync/central_america

echo 'Attempting to rsync irazu_092A_08084_081213 (255 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/irazu_092A_08084_081213 licsalert_sync/central_america

echo 'Attempting to rsync turrialba_084D_08014_091312 (256 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/turrialba_084D_08014_091312 licsalert_sync/central_america

echo 'Attempting to rsync turrialba_092A_08084_081213 (257 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/central_america/turrialba_092A_08084_081213 licsalert_sync/central_america

echo 'Attempting to rsync nevado_del_ruiz_048A_08563_131313 (258 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/nevado_del_ruiz_048A_08563_131313 licsalert_sync/south_america

echo 'Attempting to rsync nevado_del_ruiz_069D_08516_151414 (259 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/nevado_del_ruiz_069D_08516_151414 licsalert_sync/south_america

echo 'Attempting to rsync nevado_del_ruiz_150A_08574_131313 (260 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/nevado_del_ruiz_150A_08574_131313 licsalert_sync/south_america

echo 'Attempting to rsync nevado_del_huila_142D_08746_141313 (261 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/nevado_del_huila_142D_08746_141313 licsalert_sync/south_america

echo 'Attempting to rsync galeras_121A_08800_131419 (262 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/galeras_121A_08800_131419 licsalert_sync/south_america

echo 'Attempting to rsync galeras_121A_08875_090909 (263 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/galeras_121A_08875_090909 licsalert_sync/south_america

echo 'Attempting to rsync galeras_121A_08908_131313 (264 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/galeras_121A_08908_131313 licsalert_sync/south_america

echo 'Attempting to rsync galeras_121A_90526_090807 (265 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/galeras_121A_90526_090807 licsalert_sync/south_america

echo 'Attempting to rsync galeras_142D_08947_131313 (266 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/galeras_142D_08947_131313 licsalert_sync/south_america

echo 'Attempting to rsync reventador_120A_09002_080809 (267 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/reventador_120A_09002_080809 licsalert_sync/south_america

echo 'Attempting to rsync reventador_142D_08947_131313 (268 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/reventador_142D_08947_131313 licsalert_sync/south_america

echo 'Attempting to rsync guagua_pichincha_018A_09021_121620 (269 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/guagua_pichincha_018A_09021_121620 licsalert_sync/south_america

echo 'Attempting to rsync guagua_pichincha_018A_09070_111111 (270 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/guagua_pichincha_018A_09070_111111 licsalert_sync/south_america

echo 'Attempting to rsync guagua_pichincha_040D_08880_191607 (271 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/guagua_pichincha_040D_08880_191607 licsalert_sync/south_america

echo 'Attempting to rsync guagua_pichincha_040D_08919_131307 (272 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/guagua_pichincha_040D_08919_131307 licsalert_sync/south_america

echo 'Attempting to rsync guagua_pichincha_040D_09102_131313 (273 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/guagua_pichincha_040D_09102_131313 licsalert_sync/south_america

echo 'Attempting to rsync guagua_pichincha_120A_08993_090910 (274 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/guagua_pichincha_120A_08993_090910 licsalert_sync/south_america

echo 'Attempting to rsync guagua_pichincha_120A_09002_080809 (275 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/guagua_pichincha_120A_09002_080809 licsalert_sync/south_america

echo 'Attempting to rsync guagua_pichincha_121A_08908_131313 (276 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/guagua_pichincha_121A_08908_131313 licsalert_sync/south_america

echo 'Attempting to rsync guagua_pichincha_121A_09015_121211 (277 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/guagua_pichincha_121A_09015_121211 licsalert_sync/south_america

echo 'Attempting to rsync guagua_pichincha_121A_90526_090807 (278 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/guagua_pichincha_121A_90526_090807 licsalert_sync/south_america

echo 'Attempting to rsync cotopaxi_018A_09021_121620 (279 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cotopaxi_018A_09021_121620 licsalert_sync/south_america

echo 'Attempting to rsync cotopaxi_018A_09070_111111 (280 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cotopaxi_018A_09070_111111 licsalert_sync/south_america

echo 'Attempting to rsync cotopaxi_120A_08993_090910 (281 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cotopaxi_120A_08993_090910 licsalert_sync/south_america

echo 'Attempting to rsync cotopaxi_120A_09002_080809 (282 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cotopaxi_120A_09002_080809 licsalert_sync/south_america

echo 'Attempting to rsync cotopaxi_121A_09015_121211 (283 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cotopaxi_121A_09015_121211 licsalert_sync/south_america

echo 'Attempting to rsync cotopaxi_142D_09147_131313 (284 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cotopaxi_142D_09147_131313 licsalert_sync/south_america

echo 'Attempting to rsync cotopaxi_142D_09148_131313 (285 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cotopaxi_142D_09148_131313 licsalert_sync/south_america

echo 'Attempting to rsync tungurahua_018A_09021_121620 (286 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/tungurahua_018A_09021_121620 licsalert_sync/south_america

echo 'Attempting to rsync tungurahua_018A_09070_111111 (287 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/tungurahua_018A_09070_111111 licsalert_sync/south_america

echo 'Attempting to rsync tungurahua_018A_09253_131313 (288 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/tungurahua_018A_09253_131313 licsalert_sync/south_america

echo 'Attempting to rsync tungurahua_142D_09147_131313 (289 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/tungurahua_142D_09147_131313 licsalert_sync/south_america

echo 'Attempting to rsync tungurahua_142D_09148_131313 (290 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/tungurahua_142D_09148_131313 licsalert_sync/south_america

echo 'Attempting to rsync sangay_018A_09253_131313 (291 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sangay_018A_09253_131313 licsalert_sync/south_america

echo 'Attempting to rsync sangay_142D_09147_131313 (292 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sangay_142D_09147_131313 licsalert_sync/south_america

echo 'Attempting to rsync sangay_142D_09148_131313 (293 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sangay_142D_09148_131313 licsalert_sync/south_america

echo 'Attempting to rsync fernandina_106A_09048_000909 (294 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_106A_09048_000909 licsalert_sync/pacific_island

echo 'Attempting to rsync fernandina_107A_09015_000505 (295 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_107A_09015_000505 licsalert_sync/pacific_island

echo 'Attempting to rsync fernandina_107A_09018_000205 (296 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_107A_09018_000205 licsalert_sync/pacific_island

echo 'Attempting to rsync fernandina_128D_09016_110500 (297 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/fernandina_128D_09016_110500 licsalert_sync/pacific_island

echo 'Attempting to rsync wolf_106A_09048_000909 (298 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_106A_09048_000909 licsalert_sync/pacific_island

echo 'Attempting to rsync wolf_107A_09015_000505 (299 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_107A_09015_000505 licsalert_sync/pacific_island

echo 'Attempting to rsync wolf_107A_09018_000205 (300 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_107A_09018_000205 licsalert_sync/pacific_island

echo 'Attempting to rsync wolf_128D_09016_110500 (301 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/wolf_128D_09016_110500 licsalert_sync/pacific_island

echo 'Attempting to rsync sierra_negra_106A_09048_000909 (302 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/sierra_negra_106A_09048_000909 licsalert_sync/pacific_island

echo 'Attempting to rsync sierra_negra_106A_09090_000404 (303 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/sierra_negra_106A_09090_000404 licsalert_sync/pacific_island

echo 'Attempting to rsync sierra_negra_128D_09016_110500 (304 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/pacific_island/sierra_negra_128D_09016_110500 licsalert_sync/pacific_island

echo 'Attempting to rsync sabancaya_025D_10464_131313 (305 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sabancaya_025D_10464_131313 licsalert_sync/south_america

echo 'Attempting to rsync sabancaya_025D_10613_090705 (306 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sabancaya_025D_10613_090705 licsalert_sync/south_america

echo 'Attempting to rsync sabancaya_047A_10580_131313 (307 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sabancaya_047A_10580_131313 licsalert_sync/south_america

echo 'Attempting to rsync sabancaya_127D_10619_131313 (308 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/sabancaya_127D_10619_131313 licsalert_sync/south_america

echo 'Attempting to rsync ubinas_047A_10580_131313 (309 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/ubinas_047A_10580_131313 licsalert_sync/south_america

echo 'Attempting to rsync ubinas_127D_10619_131313 (310 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/ubinas_127D_10619_131313 licsalert_sync/south_america

echo 'Attempting to rsync ubinas_149A_10634_131313 (311 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/ubinas_149A_10634_131313 licsalert_sync/south_america

echo 'Attempting to rsync lascar_083D_11452_131313 (312 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/lascar_083D_11452_131313 licsalert_sync/south_america

echo 'Attempting to rsync lascar_149A_11230_131313 (313 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/lascar_149A_11230_131313 licsalert_sync/south_america

echo 'Attempting to rsync lascar_149A_11313_383737 (314 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/lascar_149A_11313_383737 licsalert_sync/south_america

echo 'Attempting to rsync lascar_149A_11428_131313 (315 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/lascar_149A_11428_131313 licsalert_sync/south_america

echo 'Attempting to rsync lascar_156D_11226_131313 (316 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/lascar_156D_11226_131313 licsalert_sync/south_america

echo 'Attempting to rsync lascar_156D_11344_313232 (317 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/lascar_156D_11344_313232 licsalert_sync/south_america

echo 'Attempting to rsync lascar_156D_11424_131313 (318 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/lascar_156D_11424_131313 licsalert_sync/south_america

echo 'Attempting to rsync planchon-peteroa_018A_12472_131313 (319 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/planchon-peteroa_018A_12472_131313 licsalert_sync/south_america

echo 'Attempting to rsync planchon-peteroa_083D_12440_131313 (320 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/planchon-peteroa_083D_12440_131313 licsalert_sync/south_america

echo 'Attempting to rsync planchon-peteroa_083D_12636_131313 (321 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/planchon-peteroa_083D_12636_131313 licsalert_sync/south_america

echo 'Attempting to rsync nevados_de_chillan_018A_12668_131313 (322 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/nevados_de_chillan_018A_12668_131313 licsalert_sync/south_america

echo 'Attempting to rsync nevados_de_chillan_083D_12636_131313 (323 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/nevados_de_chillan_083D_12636_131313 licsalert_sync/south_america

echo 'Attempting to rsync nevados_de_chillan_091A_12595_081313 (324 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/nevados_de_chillan_091A_12595_081313 licsalert_sync/south_america

echo 'Attempting to rsync nevados_de_chillan_156D_12652_131301 (325 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/nevados_de_chillan_156D_12652_131301 licsalert_sync/south_america

echo 'Attempting to rsync copahue_018A_12668_131313 (326 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/copahue_018A_12668_131313 licsalert_sync/south_america

echo 'Attempting to rsync copahue_083D_12832_131313 (327 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/copahue_083D_12832_131313 licsalert_sync/south_america

echo 'Attempting to rsync villarrica_083D_12832_131313 (328 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/villarrica_083D_12832_131313 licsalert_sync/south_america

echo 'Attempting to rsync villarrica_083D_13027_131313 (329 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/villarrica_083D_13027_131313 licsalert_sync/south_america

echo 'Attempting to rsync villarrica_164A_12971_071312 (330 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/villarrica_164A_12971_071312 licsalert_sync/south_america

echo 'Attempting to rsync calbuco_083D_13027_131313 (331 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/calbuco_083D_13027_131313 licsalert_sync/south_america

echo 'Attempting to rsync calbuco_083D_13222_131313 (332 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/calbuco_083D_13222_131313 licsalert_sync/south_america

echo 'Attempting to rsync calbuco_164A_13146_131313 (333 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/calbuco_164A_13146_131313 licsalert_sync/south_america

echo 'Attempting to rsync chaiten_083D_13222_131313 (334 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/chaiten_083D_13222_131313 licsalert_sync/south_america

echo 'Attempting to rsync chaiten_164A_13341_131313 (335 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/chaiten_164A_13341_131313 licsalert_sync/south_america

echo 'Attempting to rsync cerro_hudson_010D_13610_131313 (336 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cerro_hudson_010D_13610_131313 licsalert_sync/south_america

echo 'Attempting to rsync cerro_hudson_062A_13470_131313 (337 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cerro_hudson_062A_13470_131313 licsalert_sync/south_america

echo 'Attempting to rsync cerro_hudson_062A_13661_121313 (338 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cerro_hudson_062A_13661_121313 licsalert_sync/south_america

echo 'Attempting to rsync cerro_hudson_083D_13636_131303 (339 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cerro_hudson_083D_13636_131303 licsalert_sync/south_america

echo 'Attempting to rsync cerro_hudson_135A_13535_111313 (340 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/south_america/cerro_hudson_135A_13535_111313 licsalert_sync/south_america

echo 'Attempting to rsync soufriere_hills_054D_07353_071304 (341 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/atlantic_island/soufriere_hills_054D_07353_071304 licsalert_sync/atlantic_island

echo 'Attempting to rsync soufriere_hills_164A_07264_071317 (342 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/atlantic_island/soufriere_hills_164A_07264_071317 licsalert_sync/atlantic_island

echo 'Attempting to rsync bardarbunga_009D_02423_131316 (343 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_009D_02423_131316 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_009D_02504_202119 (344 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_009D_02504_202119 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_045A_02494_171816 (345 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_045A_02494_171816 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_045A_02538_121210 (346 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_045A_02538_121210 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_111D_02489_111313 (347 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_111D_02489_111313 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_111D_02490_152021 (348 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_111D_02490_152021 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_118A_02507_211817 (349 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_118A_02507_211817 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_147A_02466_191712 (350 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_147A_02466_191712 licsalert_sync/iceland

echo 'Attempting to rsync bardarbunga_147A_02488_131211 (351 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/iceland/bardarbunga_147A_02488_131211 licsalert_sync/iceland

echo 'Attempting to rsync terceira_002A_05136_020502 (352 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/atlantic_island/terceira_002A_05136_020502 licsalert_sync/atlantic_island

echo 'Attempting to rsync terceira_082D_05128_030500 (353 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/atlantic_island/terceira_082D_05128_030500 licsalert_sync/atlantic_island

echo 'Attempting to rsync hierro_060A_00001_030604 (354 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/atlantic_island/hierro_060A_00001_030604 licsalert_sync/atlantic_island

echo 'Attempting to rsync hierro_169D_00001_020800 (355 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/atlantic_island/hierro_169D_00001_020800 licsalert_sync/atlantic_island

echo 'Attempting to rsync fogo_075A_00001_000001 (356 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/atlantic_island/fogo_075A_00001_000001 licsalert_sync/atlantic_island

echo 'Attempting to rsync fogo_140D_07494_040200 (357 of 358)'
rsync -av --exclude '*json.gz'  --exclude 'FastICA_results.pkl' --exclude 'ICASAR_results.pkl' --exclude 'pca_results.pkl' --exclude '0*_*.png' --exclude 'mask_status.png' --exclude 'epoch_images_data.pkl' --exclude 'time_course_info.pkl' mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/atlantic_island/fogo_140D_07494_040200 licsalert_sync/atlantic_island

