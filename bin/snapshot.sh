#!/bin/bash

function refreshSnapshotNumber {

	SNAPSHOT_NUMBER=$(($SNAPSHOT_NUMBER+1))

}

if [ $# -gt 0 ]; then

	export SNAPSHOT_NUMBER

	if [ -z $SNAPSHOT_NUMBER ]; then
		SNAPSHOT_NUMBER=0
	fi

    if [ -z SURF ]; then
        export SURF=$SUBJ_DIR/'surf'
    fi

	matrixFileName="matrix.txt"
	# TODO matrix should be read from background volume
	mri_info --vox2ras $SUBJ_DIR/mri/T1.mgz --o $matrixFileName


    if [ "$1" == "use_freeview" ]; then

        screenshotName="snapshot"
	    slicesFileName="slices.txt"
	    cameraPositionsFileName="cameraPositions.txt"

	    touch $slicesFileName
	    touch $cameraPositionsFileName

        case "$2" in
		    "2vols" )
                refreshSnapshotNumber
                python bnm.recon.qc.freeview.py coronal
                freeview -v $3 $4:colormap=heat:opacity=0.3 -viewport coronal -layout 1 -cmd $slicesFileName
                python bnm.recon.qc.freeview.py sagittal
                freeview -v $3 $4:colormap=heat:opacity=0.3 -viewport sagittal -layout 1 -cmd $slicesFileName
                python bnm.recon.qc.freeview.py axial
                freeview -v $3 $4:colormap=heat:opacity=0.3 -viewport axial -layout 1 -cmd $slicesFileName
                ;;
		    "3vols" )
                refreshSnapshotNumber
                python bnm.recon.qc.freeview.py coronal
                freeview -v $3 $4:colormap=heat:opacity=0.3 $5:colormap=jet:opacity=0.5 -viewport coronal -layout 1 -cmd $slicesFileName
                python bnm.recon.qc.freeview.py sagittal
                freeview -v $3 $4:colormap=heat:opacity=0.3 $5:colormap=jet:opacity=0.5 -viewport sagittal -layout 1 -cmd $slicesFileName
                python bnm.recon.qc.freeview.py axial
                freeview -v $3 $4:colormap=heat:opacity=0.3 $5:colormap=jet:opacity=0.5 -viewport axial -layout 1 -cmd $slicesFileName
                ;;
		    "surf_annot" )
                refreshSnapshotNumber
                python bnm.recon.qc.freeview.py surface_annotation
                freeview -f $3:annot=$4 -viewport 3D -cmd $cameraPositionsFileName
                ;;
		    "vol_surf" )
                refreshSnapshotNumber
                python bnm.recon.qc.freeview.py coronal
                freeview -v $3 -f $4 -viewport coronal -cmd $slicesFileName
                python bnm.recon.qc.freeview.py sagittal
                freeview -v $3 -f $4 -viewport sagittal -cmd $slicesFileName
                python bnm.recon.qc.freeview.py axial
                freeview -v $3 -f $4 -viewport axial -cmd $slicesFileName
                ;;
		    "vol_white_pial" )
		    echo $#
                if [ $# -gt 3 and $4 == "-resampled_surface_name" ]; then
                    value=".$5"
                else
                    value=""
                fi
                refreshSnapshotNumber
                python bnm.recon.qc.freeview.py coronal
                echo $SURF/lh.white$value
                freeview -v $3 -f $SURF/{lh,rh}.{white$value:edgecolor=yellow,pial$value:edgecolor=red} -viewport coronal -cmd $slicesFileName
                python bnm.recon.qc.freeview.py sagittal
                freeview -v $3 -f $SURF/{lh,rh}.{white$value:edgecolor=yellow,pial$value:edgecolor=red} -viewport sagittal -cmd $slicesFileName
                python bnm.recon.qc.freeview.py axial
                freeview -v $3 -f $SURF/{lh,rh}.{white$value:edgecolor=yellow,pial$value:edgecolor=red} -viewport axial -cmd $slicesFileName
                ;;
	    esac
	    rm -r $slicesFileName
	    rm -r $cameraPositionsFileName

    else

        #    It needs to have bnm.recon in Python Path
        #   python setup.py develop/install
        source activate bnm-recon
        refreshSnapshotNumber
        python bnm.recon.qc.snapshot.py "$@"

    fi

    rm -r $matrixFileName

else

	echo "No arguments given!"
	exit 1

fi