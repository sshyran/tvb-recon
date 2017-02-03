#!/bin/bash

function refreshSnapshotNumber {

	SNAPSHOT_NUMBER=$(($SNAPSHOT_NUMBER+1))

}

if [ $# -gt 0 ]; then

	if [ -z $SNAPSHOT_NUMBER ]; then
			export SNAPSHOT_NUMBER=0
	fi

    if [ "$1" == "use_freeview" ]; then

        matrixFileName="matrix.txt"
	    mri_info --vox2ras $SUBJ_DIR/mri/T1.mgz --o $matrixFileName

        screenshotName="snapshot"
	    slicesFileName="slices.txt"
	    cameraPositionsFileName="cameraPositions.txt"

	    touch $slicesFileName
	    touch $cameraPositionsFileName

        case "$2" in
		    "2vols" )
                refreshSnapshotNumber
                python -m bnm.recon.qc.freeview coronal
                freeview -v $3 $4:colormap=heat:opacity=0.3 -viewport coronal -layout 1 -cmd $slicesFileName
                python -m bnm.recon.qc.freeview sagittal
                freeview -v $3 $4:colormap=heat:opacity=0.3 -viewport sagittal -layout 1 -cmd $slicesFileName
                python -m bnm.recon.qc.freeview axial
                freeview -v $3 $4:colormap=heat:opacity=0.3 -viewport axial -layout 1 -cmd $slicesFileName
                ;;
		    "3vols" )
                refreshSnapshotNumber
                python -m bnm.recon.qc.freeview coronal
                freeview -v $3 $4:colormap=heat:opacity=0.3 $5:colormap=jet:opacity=0.5 -viewport coronal -layout 1 -cmd $slicesFileName
                python -m bnm.recon.qc.freeview sagittal
                freeview -v $3 $4:colormap=heat:opacity=0.3 $5:colormap=jet:opacity=0.5 -viewport sagittal -layout 1 -cmd $slicesFileName
                python -m bnm.recon.qc.freeview axial
                freeview -v $3 $4:colormap=heat:opacity=0.3 $5:colormap=jet:opacity=0.5 -viewport axial -layout 1 -cmd $slicesFileName
                ;;
		    "surf_annot" )
                refreshSnapshotNumber
                python bnm.recon.qc.freeview.py surface_annotation
                freeview -f $3:annot=$4 -viewport 3D -cmd $cameraPositionsFileName
                ;;
		    "vol_surf" )
                refreshSnapshotNumber
                python -m bnm.recon.qc.freeview coronal
                freeview -v $3 -f $4 -viewport coronal -cmd $slicesFileName
                python -m bnm.recon.qc.freeview sagittal
                freeview -v $3 -f $4 -viewport sagittal -cmd $slicesFileName
                python -m bnm.recon.qc.freeview axial
                freeview -v $3 -f $4 -viewport axial -cmd $slicesFileName
                ;;
		    "vol_white_pial" )
		        if [ -z SURF ]; then
                    export SURF=$SUBJ_DIR/'surf'
                fi

                if [ $# -gt 3 and $4 == "-resampled_surface_name" ]; then
                    value=".$5"
                else
                    value=""
                fi
                refreshSnapshotNumber
                python -m bnm.recon.qc.freeview coronal
                freeview -v $3 -f $SURF/{lh,rh}.{white$value:edgecolor=yellow,pial$value:edgecolor=red} -viewport coronal -cmd $slicesFileName
                python -m bnm.recon.qc.freeview sagittal
                freeview -v $3 -f $SURF/{lh,rh}.{white$value:edgecolor=yellow,pial$value:edgecolor=red} -viewport sagittal -cmd $slicesFileName
                python -m bnm.recon.qc.freeview axial
                freeview -v $3 -f $SURF/{lh,rh}.{white$value:edgecolor=yellow,pial$value:edgecolor=red} -viewport axial -cmd $slicesFileName
                ;;
	    esac
	    rm -r $slicesFileName
	    rm -r $cameraPositionsFileName
	    rm -r $matrixFileName

    else

        #    It needs to have tvb.recon in Python Path
        #   python setup.py develop/install
        #   source activate tvb-recon
        refreshSnapshotNumber
        python -m bnm.recon.qc.snapshot "$@"

    fi

else

	echo "No arguments given!"
	exit 1

fi