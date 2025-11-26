""" Snasphots of the bounti segmentations

"""

import os
import sys
import nibabel as nib
import pandas as pd
import nisnap
sys.path.insert(0, os.path.abspath(os.curdir))



if __name__ == "__main__":
    dir_snap_out = "folder_to_store_the_snapshots"
    dir_fetpype = "folder_with_outputs_from_fetpype"
    subj_sess_to_process = [["subject01","session01"], ["subject02", "session01"]] # list of subjects to process
    figsize = {'x': (18, 4), 'y': (18, 4), 'z': (18, 5)}
    for subj_sess in subj_sess_to_process:
        print(subj_sess)
        subject = subj_sess[0]
        ses = subj_sess[1]
        print("________________" + subject + "_________________")
        dir_subject_out = os.path.join(dir_fetpype, subject, ses, "anat")
        recon_file = os.path.join(dir_subject_out, subject+"_"+ses+"_rec-nesvor_T2w.nii.gz")
        seg_file = os.path.join(dir_subject_out, subject+"_"+ses+"_rec-nesvor_seg-bounti_dseg.nii.gz")
        if not os.path.exists(seg_file):
            print(seg_file + " not found, skip!")
        else:
            figure = os.path.join(dir_snap_out, subject+"_"+ses + "_bounti_seg.png")
            if not os.path.exists(figure):
                done = 0
                d_max = 120
                step = 30
                while (done < 1) and (d_max > 20):
                    try:
                        slices = {'x': list(range(30, d_max, step)),
                                  'y': list(range(30, d_max, step)),
                                  'z': list(range(30, d_max, step))}
                        nisnap.plot_segment(
                            seg_file,
                            bg=recon_file,
                            slices=slices,
                            figsize=figsize,
                            # opacity=0,
                            samebox=True,
                            # labels=[1],
                            contours=True,
                            savefig=figure,
                        )
                        done = 1
                    except Exception as e:
                        print(e)
                        d_max = d_max - 20
                        step = step-5
                        print("d_max is now set to ", d_max)
