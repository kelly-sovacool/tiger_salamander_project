import plotly

plotly.offline.plot(plotly.graph_objs.Figure(data=[plotly.graph_objs.Histogram(x=[len(next(Bio.SeqIO.parse(filename, 'fasta')).seq) for filename in input.snps_filtered], name="snps filtered", autobinx=True),
                                                   plotly.graph_objs.Histogram(x=[len(next(Bio.SeqIO.parse(filename, 'fasta')).seq) for filename in input.snps_raw], name="snps unfiltered", autobinx=True)],
                                             layout=plotly.graph_objs.Layout(barmode='stack'), title=('SNP sites per locus'), xaxis=dict(title='Number of SNPs'), yaxis=dict(title='Number of loci'))),
                                             filename=output.hist, auto_open=True)
avg_snps_per_locus = sum(len(next(Bio.SeqIO.parse(filename, 'fasta')).seq) for filename in input.snps_filtered) / len(input.snps_filtered)
with open(log[0], 'w') as logfile:
    logfile.write(avg_snps_per_locus)