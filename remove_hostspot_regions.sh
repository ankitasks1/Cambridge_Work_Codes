for vars in *_peaks_id.bed
do
	echo ${vars}
	wc -l ${vars}
	bedtools intersect -a ${vars} -b merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50.txt -wa > ./hotspotfree/${vars}_hsr50.bed
done
