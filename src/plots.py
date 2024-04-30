#-------------------------------------------
#Plotting
#-------------------------------------------     


ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61


class Plot:
    @staticmethod
    def path_as_footnote(fig = None,path = "", rot = "vertical"):
        """
        Add the script origin as a footnote below or to the right side of a figure
        """        
    
        if fig == None:
            fig = plt.gcf()

        if path == None:
            import os
            try:
                path = os.path.realpath(__file__)
            except NameError: 
                # if run in Jupyter Notebook __file__ is not defined 
                # and the path is only the passed string
                path = os.path.abspath('')
        
        if rot == "vertical":
            plt.figtext(0.99, 0.5, f"{path}", horizontalalignment='right', va = "center",\
                figure = fig, rotation = "vertical", fontsize = 6)
        else:
            plt.figtext(0.99, 0, f"{path}", horizontalalignment='right', va = "bottom",\
                figure = fig, fontsize = 6)

    @staticmethod
    def add_tidal_and_f_ticks(ax,f_cpd, lower_bound = 1e-6, upper_bound = 1e-3):
        freqs = np.array(
            [24 / (14 * 24), 24 / 12.4, 2 * 24 / 12.4, 4 * 24 / 12.4, f_cpd, 2 * f_cpd, 1]
        )
        freq_labels = ["14 days", "M2", "2M2", "4M2", " \nf", " \n2f", "K1"]
        
        axis = ax.get_ylim()
        for freq in freqs:
            if freq == 24 / (14 * 24): continue
            ax.vlines(freq, lower_bound, upper_bound, color="k", alpha=0.5, linestyle="-", linewidth=0.75)
        ax.set_ylim(axis)
        
        ax2 = ax.secondary_xaxis(location="bottom")
        ax2 = Plot._axstyle(ax2, ticks="in", grid=False, spine_offset=40) 
        ax2.xaxis.set_ticks([])
        ax2.xaxis.set_ticklabels([])
        ax2.minorticks_off()
        ax2.xaxis.set_ticks(freqs)
        ax2.xaxis.set_ticklabels(freq_labels)

        return ax2

        
    @staticmethod    
    def _axstyle(
        ax=None,
        fontsize=9,
        nospine=False,
        grid=True,
        ticks="off",
        ticklength=2,
        spine_offset=5,
    ):
        """
        Apply own style to axis.

        Parameters
        ----------
        ax : AxesSubplot (optional)
            Current axis will be chosen if no axis provided

        Returns
        -------
        ax : AxesSubplot
            Axis handle
        """
        # find out background color - if this is set to ayu dark, adjust some axis
        # colors
        figcolor = plt.rcParams["figure.facecolor"]
        dark = True if figcolor == "#0d1318" else False

        if ax is None:
            ax = plt.gca()

        # Remove top and right axes lines ("spines")
        spines_to_remove = ["top", "right"]
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)

        # Remove bottom and left spines as well if desired
        if nospine:
            more_spines_to_remove = ["bottom", "left"]
            for spine in more_spines_to_remove:
                ax.spines[spine].set_visible(False)

        # For remaining spines, thin out their line and change
        # the black to a slightly off-black dark grey
        almost_black = "#262626"
        # if figure background is dark, set this close to white
        if dark:
            almost_black = "#ebe6d7"

        if ticks == "off":
            # Change the labels to the off-black
            ax.tick_params(
                axis="both",
                which="major",
                labelsize=fontsize,
                colors=almost_black,
            )
            # Get rid of ticks.
            ax.xaxis.set_ticks_position("none")
            ax.yaxis.set_ticks_position("none")
        elif ticks == "in":
            # Change the labels to the off-black
            ax.tick_params(
                axis="both",
                which="major",
                labelsize=fontsize,
                colors=almost_black,
                direction="in",
                length=ticklength,
            )

        spines_to_keep = ["bottom", "left"]
        for spine in spines_to_keep:
            ax.spines[spine].set_linewidth(0.5)
            ax.spines[spine].set_color(almost_black)
            ax.spines[spine].set_position(("outward", spine_offset))

        # Change the labels to the off-black
        ax.yaxis.label.set_color(almost_black)
        ax.yaxis.label.set_size(fontsize)
        ax.yaxis.offsetText.set_fontsize(fontsize)
        ax.xaxis.label.set_color(almost_black)
        ax.xaxis.label.set_size(fontsize)
        ax.xaxis.offsetText.set_fontsize(fontsize)

        # Change the axis title to off-black
        ax.title.set_color(almost_black)
        ax.title.set_size(fontsize + 1)

        # turn grid on
        if grid:
            ax.grid(
                which="major",
                axis="both",
                color="0.5",
                linewidth=0.25,
                linestyle="-",
                alpha=0.8,
            )

        # change legend fontsize (if there is one)
        try:
            plt.setp(ax.get_legend().get_texts(), fontsize=fontsize)
        except AttributeError:
            noleg = 1


