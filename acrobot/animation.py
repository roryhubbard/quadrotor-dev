import uuid
from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


plt.rcParams['figure.figsize'] = [15, 10]
plt.rcParams['savefig.facecolor'] = 'black'
plt.rcParams['figure.facecolor'] = 'black'
plt.rcParams['figure.edgecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'black'
plt.rcParams['axes.edgecolor'] = 'white'
plt.rcParams['axes.labelcolor'] = 'white'
plt.rcParams['axes.titlecolor'] = 'white'
plt.rcParams['xtick.color'] = 'white'
plt.rcParams['ytick.color'] = 'white'
# plt.rcParams['legend.facecolor'] = 'white'
plt.rcParams['text.color'] = 'white'
plt.rcParams["figure.autolayout"] = True


# TODO: add polar plot support
# TODO: support for collection of axes objects


def make_legend(ax, bbox_to_anchor=None, loc=None):
    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),
              bbox_to_anchor=bbox_to_anchor, loc=loc)


class AxData:

    def __init__(self, x, y, label=None, scatter=False, plot_history=True,
                 freeze_final=True, color=None, linestyle='-', marker='.',
                 markersize=None):
        self.x = x
        self.y = y
        self.label = label
        self.scatter = scatter
        self.plot_history = plot_history
        self.freeze_final = freeze_final
        self.color = color
        self.marker = marker
        self.linestyle = linestyle
        self.markersize = markersize

        if label is None:
            self.label = str(uuid.uuid4())
            self.legend = False
        else:
            self.legend = True
        if self.scatter:
            self.linestyle = ''


class Animator:

    def __init__(self, speed_multiplier, ax_data, dt, fig, ax, show_legend=True):
        self.speed_multiplier = speed_multiplier
        self.ax_data = ax_data
        self.dt = dt
        self.fig = fig
        self.ax = ax
        self.artists = dict()
        self.bbox_to_anchor = dict()
        self.loc = dict()
        self.save_path = None
        self.static = False
        self.show_legend = show_legend

    def set_leg_loc(self, bbox_to_anchor, loc='center left', ax_key='*'):
        self.bbox_to_anchor[ax_key] = bbox_to_anchor
        self.loc[ax_key] = loc

    def set_ax_limits(self, xlim=None, ylim=None, ax_key='*'):
        ax = self.get_ax_from_ax_key(ax_key)
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

    def set_save_path(self, save_path):
        self.save_path = save_path

    def make_static_plot(self):
        self.static = True

    def run(self):
        self._establish_artists()
        if self.static:
            if self.save_path is not None:
                self.fig.savefig(self.save_path)
        else:
            num_frames = self._get_num_frames()
            frame_interval = int(self.dt / self.speed_multiplier * 1000)
            anim = FuncAnimation(self.fig, self._update_frame,
                init_func=self._init_animate, frames=num_frames,
                interval=frame_interval, blit=True)
            if self.save_path is not None:
                #if '.mp4' not in self.save_path:
                #    self.save_path += '.mp4'
                if '.gif' not in self.save_path:
                    self.save_path += '.gif'
                anim.save(self.save_path)
        plt.show()
        plt.close()

    def get_ax_from_ax_key(self, ax_key):
        if ax_key == '*':
            ax = self.ax
        elif len(ax_key) > 1:
            ax = self.ax[int(ax_key[0])][int(ax_key[1])]
        else:
            ax = self.ax[int(ax_key)]
        return ax

    def _get_num_frames(self):
        num_frames = 0
        if isinstance(self.ax_data, list):
            for val in self.ax_data:
                if len(val.x) > num_frames:
                    num_frames = len(val.x)
        else:
            for values in self.ax_data.values():
                for val in values:
                    if len(val.x) > num_frames:
                        num_frames = len(val.x)
        return num_frames

    def _get_ax_artist(self, ax, label, ax_data):
        return ax.plot([], [],
            label=label,
            color=ax_data.color,
            linestyle=ax_data.linestyle,
            marker=ax_data.marker,
            markersize=ax_data.markersize,
        )

    def _add_artist(self, artist_key, ax_data, ax):
        if ax_data.legend:
            label = ax_data.label
        else:
            label = None
        self.artists[artist_key], = self._get_ax_artist(
            ax, label, ax_data)

    def _establish_artists(self):
        tmp_data = []
        if isinstance(self.ax_data, list):
            self.ax_data = {'*': self.ax_data}
        for k, v in self.ax_data.items():
            ax = self.get_ax_from_ax_key(k)
            for val in v:
                if val.scatter:
                    tmp_data.append(ax.scatter(val.x, val.y))
                else:
                    if isinstance(val.x[0], tuple) \
                            or isinstance(val.x[0], list):
                        for i in range(len(val.x)):
                            tmp_data.extend(ax.plot(val.x[i], val.y[i]))
                    else:
                        tmp_data.extend(ax.plot(val.x, val.y))
                self._add_artist(k+val.label, val, ax)
            if self.show_legend:
                make_legend(ax, self.bbox_to_anchor.get(k), self.loc.get(k))
        if not self.static:
            for _ in range(len(tmp_data)):
                tmp_data.pop().remove()

    def _init_animate(self):
        artists = []
        for artist in self.artists.values():
            artist.set_data([], [])
            artists.append(artist)
        return artists

    def _update_frame(self, i):
        artists = []
        if isinstance(self.ax_data, list):
            for val in self.ax_data:
                self._set_data(val.label, val, artists, i)
        else:
            for key, values in self.ax_data.items():
                for val in values:
                    self._set_data(key+val.label, val, artists, i)
        return artists

    def _set_data(self, artist_key, ax_data, artists, i):
        if i > len(ax_data.x) - 1:
            if ax_data.freeze_final:
                i = len(ax_data.x) - 1
            else:
                return
        if ax_data.plot_history:
            self.artists[artist_key].set_data(ax_data.x[:i+1], ax_data.y[:i+1])
        else:
            self.artists[artist_key].set_data(ax_data.x[i], ax_data.y[i])
        artists.append(self.artists[artist_key])

