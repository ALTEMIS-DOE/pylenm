def __plotUpperHalf(self, *args, **kwargs):
    corr_r = args[0].corr(args[1], 'pearson')
    corr_text = f"{corr_r:2.2f}"
    ax = plt.gca()
    ax.set_axis_off()
    marker_size = abs(corr_r) * 10000
    ax.scatter([.5], [.5], marker_size, [corr_r], alpha=0.6, cmap="coolwarm",
                vmin=-1, vmax=1, transform=ax.transAxes)
    font_size = abs(corr_r) * 40 + 5
    ax.annotate(corr_text, [.5, .48,],  xycoords="axes fraction", # [.5, .48,]
                ha='center', va='center', fontsize=font_size, fontweight='bold')


def plot_data(self, well_name, analyte_name, log_transform=True, alpha=0,
            plot_inline=True, year_interval=2, x_label='Years', y_label='', save_dir='plot_data', filter=False, col=None, equals=[]):
    """Plot concentrations over time of a specified well and analyte with a smoothed curve on interpolated data points.

    Args:

        well_name (str): name of the well to be processed
        analyte_name (str): name of the analyte to be processed
        log_transform (bool, optional): choose whether or not the data should be transformed to log base 10 values. Defaults to True.
        alpha (int, optional): alue between 0 and 10 for line smoothing. Defaults to 0.
        plot_inline (bool, optional): choose whether or not to show plot inline. Defaults to True.
        year_interval (int, optional): plot by how many years to appear in the axis e.g.(1 = every year, 5 = every 5 years, ...). Defaults to 2.
        x_label (str, optional): x axis label. Defaults to 'Years'.
        y_label (str, optional): y axis label. Defaults to ''.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_data'.
        filter (bool, optional): flag to indicate filtering. Defaults to False.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].

    Returns:
        None
    """

    # Gets appropriate data (well_name and analyte_name)
    query = self.query_data(well_name, analyte_name)
    query = self.simplify_data(data=query)

    if(type(query)==int and query == 0):
        return 'No results found for {} and {}'.format(well_name, analyte_name)
    else:
        if(filter):
            filter_res = self.filter_by_column(data=self.construction_data, col=col, equals=equals)
            if('ERROR:' in str(filter_res)):
                return filter_res
            query_wells = list(query.STATION_ID.unique())
            filter_wells = list(filter_res.index.unique())
            intersect_wells = list(set(query_wells) & set(filter_wells))
            if(len(intersect_wells)<=0):
                return 'ERROR: No results for this query with the specifed filter parameters.'
            query = query[query['STATION_ID'].isin(intersect_wells)]
        x_data = query.COLLECTION_DATE
        x_data = pd.to_datetime(x_data)
        y_data = query.RESULT
        if(log_transform):
            y_data = np.log10(y_data)

        x_RR = x_data.astype(int).to_numpy()

        nu = x_data.shape[0]

        model = SuperSmoother(alpha=alpha)
        model.fit(x_RR, y_data)
        y_pred = model.predict(x_RR)
        r = model.cv_residuals()
        out = abs(r) > 2.2*np.std(r)
        out_x = x_data[out]
        out_y = y_data[out]

        plt.figure(figsize=(8,8))
        ax = plt.axes()
        years = mdates.YearLocator(year_interval)  # every year
        months = mdates.MonthLocator()  # every month
        yearsFmt = mdates.DateFormatter('%Y')

        for label in ax.get_xticklabels():
            label.set_rotation(30)
            label.set_horizontalalignment('center')

        ax = plt.gca()
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.autoscale_view()

        unit = query.RESULT_UNITS.values[0]

        ax.set_title(str(well_name) + ' - ' + analyte_name, fontweight='bold')
        ttl = ax.title
        ttl.set_position([.5, 1.05])
        if(y_label==''):    
            if(log_transform):
                ax.set_ylabel('log-Concentration (' + unit + ')')
            else:
                ax.set_ylabel('Concentration (' + unit + ')')
        else:
            ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        small_fontSize = 15
        large_fontSize = 20
        plt.rc('axes', titlesize=large_fontSize)
        plt.rc('axes', labelsize=large_fontSize)
        plt.rc('legend', fontsize=small_fontSize)
        plt.rc('xtick', labelsize=small_fontSize)
        plt.rc('ytick', labelsize=small_fontSize) 
        ax.plot(x_data, y_data, ls='', marker='o', ms=5, color='black', alpha=1)
        ax.plot(x_data, y_pred, ls='-', marker='', ms=5, lw=2, color='blue', alpha=0.5, label="Super Smoother")
        ax.plot(out_x , out_y, ls='', marker='o', ms=5, color='red', alpha=1, label="Outliers")
        ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left', borderaxespad=0.)
        props = dict(boxstyle='round', facecolor='grey', alpha=0.15)       
        ax.text(1.05, 0.85, 'Samples: {}'.format(nu), transform=ax.transAxes, 
                fontsize=small_fontSize,
                fontweight='bold',
                verticalalignment='top', 
                bbox=props)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        # plt.savefig(save_dir + '/' + str(well_name) + '-' + analyte_name +'.png', bbox_inches="tight")
        if(plot_inline):
            plt.show()
        plt.clf()
        plt.cla()
        plt.close()


def plot_all_data(self, log_transform=True, alpha=0, year_interval=2, plot_inline=True, save_dir='plot_data'):
    """Plot concentrations over time for every well and analyte with a smoothed curve on interpolated data points.

    Args:
        log_transform (bool, optional): choose whether or not the data should be transformed to log base 10 values. Defaults to True.
        alpha (int, optional): alue between 0 and 10 for line smoothing. Defaults to 0.
        plot_inline (bool, optional): choose whether or not to show plot inline. Defaults to True.
        year_interval (int, optional): plot by how many years to appear in the axis e.g.(1 = every year, 5 = every 5 years, ...). Defaults to 2.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_data'.
    """
    analytes = ['TRITIUM','URANIUM-238','IODINE-129','SPECIFIC CONDUCTANCE', 'PH', 'DEPTH_TO_WATER']
    wells = np.array(self.data.STATION_ID.values)
    wells = np.unique(wells)
    success = 0
    errors = 0
    for well in wells:
        for analyte in analytes:
            plot = self.plot_data(well, analyte, 
                                log_transform=log_transform, 
                                alpha=alpha, 
                                year_interval=year_interval,
                                plot_inline=plot_inline,
                                save_dir=save_dir)
            if 'ERROR:' in str(plot):
                errors = errors + 1
            else:
                success = success + 1
    print("Success: ", success)
    print("Errors: ", errors)

        
def plot_MCL(self, well_name, analyte_name, year_interval=5, save_dir='plot_MCL'):
    """Plots the linear regression line of data given the analyte_name and well_name. The plot includes the prediction where the line of best fit intersects with the Maximum Concentration Limit (MCL).

    Args:
        well_name (str): ame of the well to be processed
        analyte_name (str): name of the analyte to be processed
        year_interval (int, optional): lot by how many years to appear in the axis e.g.(1 = every year, 5 = every 5 years, ...). Defaults to 5.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_MCL'.
    """
    data = self.data
    # finds the intersection point of 2 lines given the slopes and y-intercepts
    def line_intersect(m1, b1, m2, b2):
        if m1 == m2:
            print ('The lines are parallel')
            return None
        x = (b2 - b1) / (m1 - m2)
        y = m1 * x + b1
        return x,y

    # Gets appropriate data (well_name and analyte_name)
    query = self.query_data(well_name, analyte_name)

    if(type(query)==int and query == 0):
        return 'No results found for {} and {}'.format(well_name, analyte_name)
    else:   

        test = query.groupby(['COLLECTION_DATE']).mean()
        test.index = pd.to_datetime(test.index)

        x = date2num(test.index)
        y = np.log10(test.RESULT)
        ylabel = 'log-Concentration (' + self.get_unit(analyte_name) + ')'
        y = y.rename(ylabel)

        p, cov = np.polyfit(x, y, 1, cov=True)  # parameters and covariance from of the fit of 1-D polynom.

        m_unc = np.sqrt(cov[0][0])
        b_unc = np.sqrt(cov[1][1])

        f = np.poly1d(p)

        try:
            MCL = self.get_MCL(analyte_name)
            m1, b1 = f # line of best fit
            m2, b2 = 0, MCL # MCL constant

            intersection = line_intersect(m1, b1, m2, b2)

            ## Get confidence interval intersection points with MCL
            data = list(zip(x,y))
            n = len(data)
            list_slopes = []
            list_intercepts = []
            random.seed(50)
            for _ in range(80):
                sampled_data = [ random.choice(data) for _ in range(n) ]
                x_s, y_s = zip(*sampled_data)
                x_s = np.array(x_s)
                y_s = np.array(y_s)

                m_s, b_s, r, p, err = scipy.stats.linregress(x_s,y_s)
                ymodel = m_s*x_s + b_s
                list_slopes.append(m_s)
                list_intercepts.append(b_s)

            max_index = list_slopes.index(max(list_slopes))
            min_index = list_slopes.index(min(list_slopes))
            intersection_left = line_intersect(list_slopes[min_index], list_intercepts[min_index], m2, b2)
            intersection_right = line_intersect(list_slopes[max_index], list_intercepts[max_index], m2, b2)
            ##

            fig, ax = plt.subplots(figsize=(10, 6))

            ax.set_title(well_name + ' - ' + analyte_name, fontweight='bold')
            ttl = ax.title
            ttl.set_position([.5, 1.05])
            years = mdates.YearLocator(year_interval)  # every year
            months = mdates.MonthLocator()  # every month
            yearsFmt = mdates.DateFormatter('%Y') 

            ax.xaxis.set_major_locator(years)
            ax = plt.gca()
            ax.xaxis.set_major_locator(years)
            ax.xaxis.set_major_formatter(yearsFmt)
            ax.autoscale_view()
            ax.grid(True, alpha=0.4)
            small_fontSize = 15
            large_fontSize = 20
            plt.rc('axes', titlesize=large_fontSize)
            plt.rc('axes', labelsize=large_fontSize)
            plt.rc('legend', fontsize=small_fontSize)
            plt.rc('xtick', labelsize=small_fontSize)
            plt.rc('ytick', labelsize=small_fontSize)

            ax.set_xlabel('Years')
            ax.set_ylabel('log-Concentration (' + self.get_unit(analyte_name) + ')')

            if(intersection[0] < min(x)):
                temp = intersection_left
                intersection_left = intersection_right
                intersection_right = temp
                ax.set_ylim([0, max(y)+1])
                ax.set_xlim([intersection_left[0]-1000, max(x)+1000])
            elif(intersection[0] < max(x) and intersection[0] > min(x)):
                ax.set_ylim([0, max(y)+1])
                ax.set_xlim(min(x)-1000, max(x)+1000)

            else:
                ax.set_ylim([0, max(y)+1])
                ax.set_xlim([min(x)-1000, intersection_right[0]+1000])

            ax = sns.regplot(x, y, logx=True, truncate=False, seed=42, n_boot=1000, ci=95) # Line of best fit
            ax.plot(x, y, ls='', marker='o', ms=5, color='black', alpha=1) # Data
            ax.axhline(y=MCL, color='r', linestyle='--') # MCL
            ax.plot(intersection[0], intersection[1], color='blue', marker='o', ms=10)
            ax.plot(intersection_left[0], intersection_left[1], color='green', marker='o', ms=5)
            ax.plot(intersection_right[0], intersection_right[1], color='green', marker='o', ms=5)

            predict = num2date(intersection[0]).date()
            l_predict = num2date(intersection_left[0]).date()
            u_predict = num2date(intersection_right[0]).date()
            ax.annotate(predict, (intersection[0], intersection[1]), xytext=(intersection[0], intersection[1]+1), 
                        bbox=dict(boxstyle="round", alpha=0.1),ha='center', arrowprops=dict(arrowstyle="->", color='blue'), fontsize=small_fontSize, fontweight='bold')
            props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
            ax.text(1.1, 0.5, 'Lower confidence:  {}\n            Prediction:  {}\nUpper confidence:  {}'.format(l_predict, predict, u_predict), transform=ax.transAxes, fontsize=small_fontSize, fontweight='bold', verticalalignment='bottom', bbox=props)

            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            plt.savefig(save_dir + '/' + well_name + '-' + analyte_name +'.png', bbox_inches="tight")

        except:
            print('ERROR: Something went wrong')
            return None

