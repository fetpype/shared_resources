""" split the segmentation map in two hemispheres
"""
import sys; print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['../fet-processing'])

import os
import pickle
import ants

from configuration import (
    MarsFet_HASTE_NORMATIVE_SUB_SESS_PICKL,
    MarsFet_OUPUT_NNUNET_SEGMENTATION_PATH,
    MarsFet_SRR_NESVOR_PATH,
    ATLAS_GHOLIPOUR_PATH
)


def ants_register(fixed, atlas_week):
    moving_atlas_file = os.path.join(ATLAS_GHOLIPOUR_PATH, "STA" + atlas_week + "_flipped.nii.gz")
    # load atlas volume
    moving_atlas = ants.image_read(moving_atlas_file)
    # fixed.plot(overlay=moving_atlas, title='Before Registration', overlay_alpha = 0.5)
    # comute registration
    mytx = ants.registration(fixed=fixed, moving=moving_atlas, type_of_transform="Affine")  # 'SyN' )
    # fwdtransforms: Transforms to move from moving to fixed image.
    # invtransforms: Transforms to move from fixed to moving image.
    # fwdtransform = mytx['fwdtransforms']
    warped_atlas = mytx['warpedmovout']
    # warped_atlas.plot()
    # compute MI to find the closest atlas
    wraped_mi = ants.image_mutual_information(fixed, warped_atlas)
    return wraped_mi


if __name__ == "__main__":
    sequence = "haste"
    with open(MarsFet_HASTE_NORMATIVE_SUB_SESS_PICKL, 'rb') as file:
        subj_sess = pickle.load(file)
    print(str(len(subj_sess))+' subjects to be processed')
    dir_segmentation = MarsFet_OUPUT_NNUNET_SEGMENTATION_PATH

    failed_subj = list()
    for subj_s in subj_sess:
        subj_s_split = subj_s.split("/")
        subject = subj_s_split[0]
        session = subj_s_split[1]
        print(subject+'_'+session)
        file_seg_out = os.path.join(dir_segmentation, subject + "_" + session + "_" + "acq-" + sequence + "_rec-nesvor_desc-aligned_T2w_seg_hemi.nii.gz")
        #file_seg_out = "/home/toz/test_hemi_split/output/sub-0001_ses-0001_acq-haste_rec-nesvor_desc-aligned_T2w_split.nii.gz"
        if os.path.exists(file_seg_out):
            print('already processed, skip')
        else:
            file_seg_in = os.path.join(dir_segmentation, subject + "_" + session + "_" + "acq-" + sequence + "_rec-nesvor_desc-aligned_T2w.nii.gz")
            #file_seg_in = "/home/toz/test_hemi_split/sub-0001/ses-0001/sub-0001_ses-0001_acq-haste_rec-nesvor_desc-aligned_T2w.nii.gz"

            dir_reconst = os.path.join(
                MarsFet_SRR_NESVOR_PATH,
                subj_s,
                sequence,
            )

            file_T2_subj = os.path.join(
                dir_reconst,
                "default_reconst",
                subject + "_" + session + "_" + "acq-" + sequence + "_rec-nesvor_desc-aligned_T2w.nii.gz",
            )
            #file_T2_subj = "/home/toz/test_hemi_split/sub-0001/ses-0001/haste/default_reconst/sub-0001_ses-0001_acq-haste_rec-nesvor_desc-aligned_T2w.nii.gz"

            try:
                # load individual volumes
                fixed = ants.image_read(file_T2_subj)
                fixed_seg = ants.image_read(file_seg_in)
                # find the closest atlas
                atlas_mi = list()
                wraped_mi_29 = ants_register(fixed, "29")
                print(["29", wraped_mi_29])
                wraped_mi_30 = ants_register(fixed, "30")
                print(["30", wraped_mi_30])
                diff_wraped_mi = wraped_mi_30 - wraped_mi_29
                print("30 - 29 = "+str(diff_wraped_mi))
                if diff_wraped_mi<0: # atlas 30 closer than atlas 29 so test older atlases
                    atlas_mi.append(["29", wraped_mi_29])
                    atlas_mi.append(["30", wraped_mi_30])
                    best_atlas = ["30", wraped_mi_30]
                    atlas_week_num = 30
                    while diff_wraped_mi<0 and atlas_week_num<38:
                        atlas_week_num += 1
                        atlas_week = str(atlas_week_num)
                        wraped_mi = ants_register(fixed, atlas_week)
                        atlas_mi.append([atlas_week, wraped_mi])
                        print(atlas_mi[-1])
                        if wraped_mi < best_atlas[1]:
                            best_atlas = [atlas_week, wraped_mi]
                        diff_wraped_mi = atlas_mi[-1][1] - atlas_mi[-2][1]
                        print(atlas_mi[-1][0]+" - "+atlas_mi[-2][0]+" = "+str(diff_wraped_mi))
                else: # atlas 30 closer than atlas 29 so test younger atlases
                    atlas_mi.append(["30", wraped_mi_30])
                    atlas_mi.append(["29", wraped_mi_29])
                    best_atlas = ["29", wraped_mi_29]
                    diff_wraped_mi = -diff_wraped_mi
                    atlas_week_num = 29
                    while diff_wraped_mi<0 and atlas_week_num>21:
                        atlas_week_num -= 1
                        atlas_week = str(atlas_week_num)
                        wraped_mi = ants_register(fixed, atlas_week)
                        atlas_mi.append([atlas_week, wraped_mi])
                        print(atlas_mi[-1])
                        if wraped_mi < best_atlas[1]:
                            best_atlas = [atlas_week, wraped_mi]
                        diff_wraped_mi = atlas_mi[-1][1] - atlas_mi[-2][1]
                        print(atlas_mi[-1][0]+" - "+atlas_mi[-2][0]+" = "+str(diff_wraped_mi))
                print("______________________")
                print("The closest atlas is "+best_atlas[0])
                print("Compute non-linear registration")
                # non-linear registration on the closest atlas
                best_atlas_week = best_atlas[0]
                best_atlas_file = os.path.join(ATLAS_GHOLIPOUR_PATH, "STA" + best_atlas_week + "_flipped.nii.gz")
                # load best atlas volume and corresponding hemi seg
                moving_best_atlas = ants.image_read(best_atlas_file)
                moving_seg_file = os.path.join(ATLAS_GHOLIPOUR_PATH, "STA" + best_atlas_week + "_all_reg_LR_dilall_flipped.nii.gz")
                moving_best_seg = ants.image_read(moving_seg_file)
                # non-linear registration of the best atlas
                mytx_best = ants.registration(fixed=fixed , moving=moving_best_atlas, type_of_transform='SyN' )
                # fwdtransforms: Transforms to move from moving to fixed image.
                # invtransforms: Transforms to move from fixed to moving image.
                fwdtransform_best = mytx_best['fwdtransforms']
                warped_best_atlas = mytx_best['warpedmovout']
                #fixed.plot(overlay=warped_atlas,
                #           title='After Registration', overlay_alpha = 0.5)
                wraped_mi = ants.image_mutual_information( fixed, warped_best_atlas )
                atlas_mi.append([best_atlas_week,wraped_mi])
                print(atlas_mi[-1])

                warped_best_seg = ants.apply_transforms(fixed=fixed, moving=moving_best_seg,
                                                transformlist=mytx_best['fwdtransforms'],
                                                interpolator="nearestNeighbor")
                #ants.image_write(warped_best_atlas, aligned_image_file)
                #warpedimage.plot()
                # fixed.plot(overlay=warped_seg, title='seg on fixed', overlay_alpha=0.5)

                ## use the aligned atlas hemi to split the segmentation
                seg_arr = fixed_seg.numpy()
                bm = seg_arr>0
                seg_atlas = warped_best_seg.numpy()
                out_arr = seg_arr + 100*seg_atlas
                out_arr = out_arr*bm
                # recombine brainstem
                brainstem = seg_arr == 8
                out_arr[brainstem] = 8
                # recombine background
                brainstem = seg_arr == 4
                out_arr[brainstem] = 4
                seg_out = fixed_seg.new_image_like(out_arr)
                ants.image_write(seg_out, file_seg_out)
                print("splitted segmentation saved as:")
                print(file_seg_out)
            except Exception as e:
                print(e)
                failed_subj.append(subj_s)
    print(str(len(failed_subj))+' failed subjects')
    print(failed_subj)


