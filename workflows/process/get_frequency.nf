process get_frequency {
        // errorStrategy is set in the nodes.config file
        label 'frequency'
        publishDir "${params.output}/${name}/pileup/", mode: 'copy'
    input:
        tuple val(name), path(fasta), path(positions), path(bam) 

  	output:
    	tuple val(name), file("${name}_Frequency.tsv"), path("${name}_pileup.tsv")
  	script:
    if (task.exitStatus == 137) 

    """
    # ACTIVATE HISTORY
        set -euxo pipefail

    # remove headers and create one correct bed position file for samtools (name \t position)
    tail -q -n+2 ${positions} | awk -v OFS='\\t' '{print \$1,\$5}' > clean.bed

    # pileup for all positions
    samtools index ${bam}
    samtools mpileup ${bam} --max-depth 200 --no-BAQ -f ${fasta} --positions clean.bed --output ${name}_pileup.tsv


    # count stuff per degenerate base
        awk -F'\\t' '{
            # Process Y
            if (\$3 == "Y") {
                total_Y_upper += gsub(/[CT]/, "&", \$5);
                count_YT_upper += gsub(/T/, "&", \$5);
                count_YC_upper += gsub(/C/, "&", \$5);
                total_y_lower += gsub(/[ct]/, "&", \$5);
                count_Yt_lower += gsub(/t/, "&", \$5);
                count_Yc_lower += gsub(/c/, "&", \$5);
            }
            # Process R
            if (\$3 == "R") {
                total_R_upper += gsub(/[AG]/, "&", \$5);
                count_RA_upper += gsub(/A/, "&", \$5);
                count_RG_upper += gsub(/G/, "&", \$5);
                total_r_lower += gsub(/[ag]/, "&", \$5);
                count_Ra_lower += gsub(/a/, "&", \$5);
                count_Rg_lower += gsub(/g/, "&", \$5);
            }
            # Process W
            if (\$3 == "W") {
                total_W_upper += gsub(/[AT]/, "&", \$5);
                count_WA_upper += gsub(/A/, "&", \$5);
                count_WT_upper += gsub(/T/, "&", \$5);
                total_w_lower += gsub(/[at]/, "&", \$5);
                count_Wa_lower += gsub(/a/, "&", \$5);
                count_Wt_lower += gsub(/t/, "&", \$5);
            }
            # Process S
            if (\$3 == "S") {
                total_S_upper += gsub(/[GC]/, "&", \$5);
                count_SG_upper += gsub(/G/, "&", \$5);
                count_SC_upper += gsub(/C/, "&", \$5);
                total_s_lower += gsub(/[gc]/, "&", \$5);
                count_Sg_lower += gsub(/g/, "&", \$5);
                count_Sc_lower += gsub(/c/, "&", \$5);
            }
            # Process M
            if (\$3 == "M") {
                total_M_upper += gsub(/[AC]/, "&", \$5);
                count_MA_upper += gsub(/A/, "&", \$5);
                count_MC_upper += gsub(/C/, "&", \$5);
                total_m_lower += gsub(/[ac]/, "&", \$5);
                count_Ma_lower += gsub(/a/, "&", \$5);
                count_Mc_lower += gsub(/c/, "&", \$5);
            }
            # Process K
            if (\$3 == "K") {
                total_K_upper += gsub(/[TG]/, "&", \$5);
                count_KT_upper += gsub(/T/, "&", \$5);
                count_KG_upper += gsub(/G/, "&", \$5);
                total_k_lower += gsub(/[tg]/, "&", \$5);
                count_Kt_lower += gsub(/t/, "&", \$5);
                count_Kg_lower += gsub(/g/, "&", \$5);
            }
        }
        END {
            # Print results for Y
            if (total_Y_upper > 0) print "forward\\tY\\tT\\t" count_YT_upper / total_Y_upper >> "${name}_Frequency.tsv";
            if (total_Y_upper > 0) print "forward\\tY\\tC\\t" count_YC_upper / total_Y_upper >> "${name}_Frequency.tsv";
            if (total_y_lower > 0) print "reverse\\tY\\tT\\t" count_Yt_lower / total_y_lower >> "${name}_Frequency.tsv";
            if (total_y_lower > 0) print "reverse\\tY\\tC\\t" count_Yc_lower / total_y_lower >> "${name}_Frequency.tsv";

            # Print results for R
            if (total_R_upper > 0) print "forward\\tR\\tA\\t" count_RA_upper / total_R_upper >> "${name}_Frequency.tsv";
            if (total_R_upper > 0) print "forward\\tR\\tG\\t" count_RG_upper / total_R_upper >> "${name}_Frequency.tsv";
            if (total_r_lower > 0) print "reverse\\tR\\tA\\t" count_Ra_lower / total_r_lower >> "${name}_Frequency.tsv";
            if (total_r_lower > 0) print "reverse\\tR\\tG\\t" count_Rg_lower / total_r_lower >> "${name}_Frequency.tsv";

            # Print results for W
            if (total_W_upper > 0) print "forward\\tW\\tA\\t" count_WA_upper / total_W_upper >> "${name}_Frequency.tsv";
            if (total_W_upper > 0) print "forward\\tW\\tT\\t" count_WT_upper / total_W_upper >> "${name}_Frequency.tsv";
            if (total_w_lower > 0) print "reverse\\tW\\tA\\t" count_Wa_lower / total_w_lower >> "${name}_Frequency.tsv";
            if (total_w_lower > 0) print "reverse\\tW\\tT\\t" count_Wt_lower / total_w_lower >> "${name}_Frequency.tsv";

            # Print results for S
            if (total_S_upper > 0) print "forward\\tS\\tG\\t" count_SG_upper / total_S_upper >> "${name}_Frequency.tsv";
            if (total_S_upper > 0) print "forward\\tS\\tC\\t" count_SC_upper / total_S_upper >> "${name}_Frequency.tsv";
            if (total_s_lower > 0) print "reverse\\tS\\tG\\t" count_Sg_lower / total_s_lower >> "${name}_Frequency.tsv";
            if (total_s_lower > 0) print "reverse\\tS\\tC\\t" count_Sc_lower / total_s_lower >> "${name}_Frequency.tsv";

            # Print results for M
            if (total_M_upper > 0) print "forward\\tM\\tA\\t" count_MA_upper / total_M_upper >> "${name}_Frequency.tsv";
            if (total_M_upper > 0) print "forward\\tM\\tC\\t" count_MC_upper / total_M_upper >> "${name}_Frequency.tsv";
            if (total_m_lower > 0) print "reverse\\tM\\tA\\t" count_Ma_lower / total_m_lower >> "${name}_Frequency.tsv";
            if (total_m_lower > 0) print "reverse\\tM\\tC\\t" count_Mc_lower / total_m_lower >> "${name}_Frequency.tsv";

            # Print results for K
            if (total_K_upper > 0) print "forward\\tK\\tG\\t" count_KG_upper / total_K_upper >> "${name}_Frequency.tsv";
            if (total_K_upper > 0) print "forward\\tK\\tT\\t" count_KT_upper / total_K_upper >> "${name}_Frequency.tsv";
            if (total_k_lower > 0) print "reverse\\tK\\tG\\t" count_Kg_lower / total_k_lower >> "${name}_Frequency.tsv";
            if (total_k_lower > 0) print "reverse\\tK\\tT\\t" count_Kt_lower / total_k_lower >> "${name}_Frequency.tsv";
        }' "${name}_pileup.tsv"

    """
    else
    """
        # ACTIVATE HISTORY
        set -euxo pipefail

    # remove headers and create one correct bed position file for samtools (name \t position)
    tail -q -n+2 ${positions} | awk -v OFS='\\t' '{print \$1,\$5}' > clean.bed

    # pileup for all positions
    samtools index ${bam}
    samtools mpileup ${bam} --max-depth 200 -f ${fasta} --positions clean.bed --output ${name}_pileup.tsv


    # count stuff per degenerate base
        awk -F'\\t' '{
            # Process Y
            if (\$3 == "Y") {
                total_Y_upper += gsub(/[CT]/, "&", \$5);
                count_YT_upper += gsub(/T/, "&", \$5);
                count_YC_upper += gsub(/C/, "&", \$5);
                total_y_lower += gsub(/[ct]/, "&", \$5);
                count_Yt_lower += gsub(/t/, "&", \$5);
                count_Yc_lower += gsub(/c/, "&", \$5);
            }
            # Process R
            if (\$3 == "R") {
                total_R_upper += gsub(/[AG]/, "&", \$5);
                count_RA_upper += gsub(/A/, "&", \$5);
                count_RG_upper += gsub(/G/, "&", \$5);
                total_r_lower += gsub(/[ag]/, "&", \$5);
                count_Ra_lower += gsub(/a/, "&", \$5);
                count_Rg_lower += gsub(/g/, "&", \$5);
            }
            # Process W
            if (\$3 == "W") {
                total_W_upper += gsub(/[AT]/, "&", \$5);
                count_WA_upper += gsub(/A/, "&", \$5);
                count_WT_upper += gsub(/T/, "&", \$5);
                total_w_lower += gsub(/[at]/, "&", \$5);
                count_Wa_lower += gsub(/a/, "&", \$5);
                count_Wt_lower += gsub(/t/, "&", \$5);
            }
            # Process S
            if (\$3 == "S") {
                total_S_upper += gsub(/[GC]/, "&", \$5);
                count_SG_upper += gsub(/G/, "&", \$5);
                count_SC_upper += gsub(/C/, "&", \$5);
                total_s_lower += gsub(/[gc]/, "&", \$5);
                count_Sg_lower += gsub(/g/, "&", \$5);
                count_Sc_lower += gsub(/c/, "&", \$5);
            }
            # Process M
            if (\$3 == "M") {
                total_M_upper += gsub(/[AC]/, "&", \$5);
                count_MA_upper += gsub(/A/, "&", \$5);
                count_MC_upper += gsub(/C/, "&", \$5);
                total_m_lower += gsub(/[ac]/, "&", \$5);
                count_Ma_lower += gsub(/a/, "&", \$5);
                count_Mc_lower += gsub(/c/, "&", \$5);
            }
            # Process K
            if (\$3 == "K") {
                total_K_upper += gsub(/[TG]/, "&", \$5);
                count_KT_upper += gsub(/T/, "&", \$5);
                count_KG_upper += gsub(/G/, "&", \$5);
                total_k_lower += gsub(/[tg]/, "&", \$5);
                count_Kt_lower += gsub(/t/, "&", \$5);
                count_Kg_lower += gsub(/g/, "&", \$5);
            }
        }
        END {
            # Print results for Y
            if (total_Y_upper > 0) print "forward\\tY\\tT\\t" count_YT_upper / total_Y_upper >> "${name}_Frequency.tsv";
            if (total_Y_upper > 0) print "forward\\tY\\tC\\t" count_YC_upper / total_Y_upper >> "${name}_Frequency.tsv";
            if (total_y_lower > 0) print "reverse\\tY\\tT\\t" count_Yt_lower / total_y_lower >> "${name}_Frequency.tsv";
            if (total_y_lower > 0) print "reverse\\tY\\tC\\t" count_Yc_lower / total_y_lower >> "${name}_Frequency.tsv";

            # Print results for R
            if (total_R_upper > 0) print "forward\\tR\\tA\\t" count_RA_upper / total_R_upper >> "${name}_Frequency.tsv";
            if (total_R_upper > 0) print "forward\\tR\\tG\\t" count_RG_upper / total_R_upper >> "${name}_Frequency.tsv";
            if (total_r_lower > 0) print "reverse\\tR\\tA\\t" count_Ra_lower / total_r_lower >> "${name}_Frequency.tsv";
            if (total_r_lower > 0) print "reverse\\tR\\tG\\t" count_Rg_lower / total_r_lower >> "${name}_Frequency.tsv";

            # Print results for W
            if (total_W_upper > 0) print "forward\\tW\\tA\\t" count_WA_upper / total_W_upper >> "${name}_Frequency.tsv";
            if (total_W_upper > 0) print "forward\\tW\\tT\\t" count_WT_upper / total_W_upper >> "${name}_Frequency.tsv";
            if (total_w_lower > 0) print "reverse\\tW\\tA\\t" count_Wa_lower / total_w_lower >> "${name}_Frequency.tsv";
            if (total_w_lower > 0) print "reverse\\tW\\tT\\t" count_Wt_lower / total_w_lower >> "${name}_Frequency.tsv";

            # Print results for S
            if (total_S_upper > 0) print "forward\\tS\\tG\\t" count_SG_upper / total_S_upper >> "${name}_Frequency.tsv";
            if (total_S_upper > 0) print "forward\\tS\\tC\\t" count_SC_upper / total_S_upper >> "${name}_Frequency.tsv";
            if (total_s_lower > 0) print "reverse\\tS\\tG\\t" count_Sg_lower / total_s_lower >> "${name}_Frequency.tsv";
            if (total_s_lower > 0) print "reverse\\tS\\tC\\t" count_Sc_lower / total_s_lower >> "${name}_Frequency.tsv";

            # Print results for M
            if (total_M_upper > 0) print "forward\\tM\\tA\\t" count_MA_upper / total_M_upper >> "${name}_Frequency.tsv";
            if (total_M_upper > 0) print "forward\\tM\\tC\\t" count_MC_upper / total_M_upper >> "${name}_Frequency.tsv";
            if (total_m_lower > 0) print "reverse\\tM\\tA\\t" count_Ma_lower / total_m_lower >> "${name}_Frequency.tsv";
            if (total_m_lower > 0) print "reverse\\tM\\tC\\t" count_Mc_lower / total_m_lower >> "${name}_Frequency.tsv";

            # Print results for K
            if (total_K_upper > 0) print "forward\\tK\\tG\\t" count_KG_upper / total_K_upper >> "${name}_Frequency.tsv";
            if (total_K_upper > 0) print "forward\\tK\\tT\\t" count_KT_upper / total_K_upper >> "${name}_Frequency.tsv";
            if (total_k_lower > 0) print "reverse\\tK\\tG\\t" count_Kg_lower / total_k_lower >> "${name}_Frequency.tsv";
            if (total_k_lower > 0) print "reverse\\tK\\tT\\t" count_Kt_lower / total_k_lower >> "${name}_Frequency.tsv";
        }' "${name}_pileup.tsv"
    """
}

