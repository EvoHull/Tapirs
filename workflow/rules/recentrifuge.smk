# =====================================================
# KRAKEN2 RECENTRIFUGE
# =====================================================

# -----------------------------------------------------
# Recentrifuge produces interactive html displays
# of kraken output
# -----------------------------------------------------

rule kraken2_recentrifuge:
    conda:
        "envs/tax.yaml"
    input:
        taxdb = config["taxdump"],
        kraken_output = expand("results/kraken2/outputs/{real_combos}.krk", real_combos = real_combos)
    params:
        directory("results/kraken2/outputs/")
    output:
        re_html = "results/recentrifuge/" + config['my_experiment'] + ".html",
        re_xlsx = "results/recentrifuge/" + config['my_experiment'] + ".xlsx"
    shell:
        "rcf -a \
        -n {input.taxdb} \
        -k {params} \
        -o {output.re_html}"



# rule kraken2_recentrifuge:
#     conda:
#         config['conda']
#     input:
#         taxdb = config["taxdump"],
#         kraken_out_N = "results/kraken/{LIBRARIES}/{SAMPLES}.krk",
#         rcf = "recentrifuge/retaxdump"
#     output:
#         "reports/recentrifuge/{LIBRARIES}/{SAMPLES}.krk.html"
#     shell:
#         "./recentrifuge/rcf -n {input.taxdb} -k {input.kraken_out_N} -o {output}"
