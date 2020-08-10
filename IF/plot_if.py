import numpy as np
import bokeh
from bokeh.models import ColumnDataSource
from utils import Surface3d
from compute_impact import ImpactFactor
from copy import deepcopy


def plot_chunks(out_dir):
    for bs_num in np.arange(1, 6):
        for bs_idx in range(bs_num):
            for prev_change in np.arange(0, bs_idx+1, step=0.1):
                if prev_change > bs_idx:
                    pass
                bs_change = np.zeros(bs_num)
                new_idx = int(np.floor(prev_change))
                bs_change[:new_idx] = 1
                bs_change[new_idx] = (prev_change) % 1

                for strategy in ImpactFactor.available_strategies:
                    imp = ImpactFactor(strategy=strategy)

                    expression = np.arange(1, -0.05, -0.1)
                    single_bs_change = np.array([1])  # np.arange(0, 1, 0.1)
                    xx, yy = np.meshgrid(expression, single_bs_change)
                    xx = xx.ravel()
                    yy = yy.ravel()

                    tissue = []
                    transcript = []
                    gene = []
                    for x, y in zip(xx, yy):
                        changes = deepcopy(bs_change)
                        changes[bs_idx] = y
                        tissue.append(imp.transcript_in_tissiue_if(
                            x, 3, 5, changes))
                        transcript.append(imp.transcript_if(
                            x, 7, 3, 5, changes))
                        gene.append(imp.gene_if(x, 7, 3, changes))

                    values = {'tissue': tissue,
                              'transcript': transcript,
                              'gene': gene
                              }

                    bs_change[bs_idx] = -1
                    print(bs_change)
                    bs_string = '_'.join(str(bs_change))
                    for level in values:
                        out = '{}/{}-{}-{}_BS-{}.html'.format(out_dir,
                                                              strategy, level, bs_num, bs_string)

                        source = ColumnDataSource(
                            data=dict(x=xx, y=yy, z=values[level]))
                        surface = Surface3d(x="x", y="y", z="z",
                                            data_source=source, width=1600, height=1600)
                        bokeh.io.output_file(out)
                        bokeh.io.save(surface)


def plot_landscapes(out_dir):
    for strategy in ImpactFactor.available_strategies:
        imp = ImpactFactor(strategy=strategy)
        for bs_num in np.arange(1, 6):
            landscape_x = []
            landscape_y = []
            values = {'tissue': [],
                      'transcript': [],
                      'gene': []
                      }
            expression = np.arange(0, 1, 0.1)
            single_bs_change = np.array([1])  # np.arange(0, 1, 0.1)
            xx, yy = np.meshgrid(expression, single_bs_change)
            x_coef = 30
            y_coef = 1
            xx = xx.ravel()*x_coef
            yy = yy.ravel()*y_coef
            bs_change = None

            for bs_idx in range(bs_num):
                for prev_change in np.arange(0, bs_idx+1, step=0.1):
                    landscape_x += list(xx)  # landscape_x + [a + (prev_change*10)*max(xx) for a in xx][::-1]
                    landscape_y = landscape_y + [a + (prev_change*10)*max(yy) for a in yy][::-1]  # list(yy)
                    bs_change = np.zeros(bs_num)
                    new_idx = int(np.floor(prev_change))
                    bs_change[:new_idx] = 1
                    bs_change[new_idx] = (prev_change) % 1
                    tissue = []
                    transcript = []
                    gene = []
                    print(bs_idx, new_idx)
                    print(bs_change)
                    for x, y in zip(xx, yy):
                        tissue.append(imp.transcript_in_tissue_if(x, 3, 5, bs_change))
                        transcript.append(imp.transcript_if(x, 7, 3, 5, bs_change))
                        gene.append(imp.gene_if(x, 7, 3, bs_change))

                    values['tissue'] += tissue
                    values['transcript'] += transcript
                    values['gene'] += gene
            print(bs_change)
            for level in values:
                values[level] = values[level]/max(values[level])
                out = '{}/{}-{}-{}_BS.html'.format(out_dir, strategy, level, bs_num)
                source = ColumnDataSource(data=dict(x=landscape_x, y=landscape_y, z=values[level]))
                surface = Surface3d(x="x", y="y", z="z", data_source=source)
                bokeh.io.output_file(out)
                bokeh.io.save(surface)


if __name__ == "__main__":
    out_dir = 'figures'
    plot_landscapes(out_dir)

# print('Tissue: {}\nTranscript: {}\n Gene: {}'.format(tissue, transcript, gene))
