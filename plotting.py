from matplotlib import pyplot as plt

def make_barplot(xlabels, values, xlabel = None, ylabel = None, title = None, width = 0.5, ax = None, logscale = False):
    if ax is None:
        fig, ax = plt.subplots()
    margin = (1 - width) + width / 2
    ax.set_xlim(-margin, len(xlabels) - 1 + margin)
    if ylabel:
        ax.set(ylabel=ylabel)
    if xlabel:
        ax.set(xlabel=xlabel)
    if logscale:
        ax.set_yscale("log")
    if title:
        ax.set_title(title)
    ax.bar(xlabels, values, width=width, align="center")
