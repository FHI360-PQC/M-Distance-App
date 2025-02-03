import plotly.graph_objects as go
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shinywidgets import output_widget, render_widget  
from shiny.types import FileInfo
from PQCDRS import *

file_formats = {"tellspec": "Tellspec", "labspec": "Labspec", 
                "ti": "Texas Instruments", "compressed": "Compressed"}
default_format = "compressed"
accepted_extensions = ['.csv','.txt']

### UI ####
app_ui = ui.page_fluid(
    ui.navset_card_pill(
        ### Reference Upload
        ui.nav_panel("Reference",
            ui.card(
                ui.layout_sidebar(
                    ui.sidebar(
                    ui.input_select("reftype", "Reference File Format:", file_formats,
                        selected=default_format
                    ),  
                    ui.input_file('refdata', 'Upload Reference Spectra',
                        multiple=False, accept=accepted_extensions
                    ),
                    ui.input_checkbox("refbkg_check", "Apply Background?", False),
                    ui.panel_conditional("input.refbkg_check",
                        ui.input_file('refbkg', 'Upload Reference Background',
                            multiple=False, accept=accepted_extensions
                        )
                    ),
                    ui.input_checkbox("ref_logtrans", "Log Transform", False),
                    ui.input_checkbox("ref_trunc", "Truncation", False),
                    ui.panel_conditional("input.ref_trunc",
                        ui.input_numeric("ref_truncmin", "Truncation Lowerbound", 0, min=0, max=1, step=1),
                        ui.input_numeric("ref_truncmax", "Truncation Upperbound", 1, min=0, max=1, step=1),
                    ),
                    bg="#f8f8f8", open='always'
                    ), # end ui.sidebar
                    ui.navset_card_pill(
                        ui.nav_panel("Spectra",
                            ui.output_plot("Ref_SpectraPlot"),
                            ui.output_data_frame("Ref_SpectraTable"),
                        ),
                        ui.nav_panel("PCA",
                            ui.input_checkbox("pca_center", "Center", True),
                            ui.input_checkbox("pca_scale", "Scale", False),
                            ui.input_numeric("pca_ncomp", "# of Components", 5, min=1, step=1),
                            ui.input_numeric("mdist_ci", "Confidence Interval", 0.95, min=0, max=1, step=0.01),
                            output_widget("PCA_ScreePlot"),
                            ui.output_plot("PCA_ResidualPlot"),
                            ui.output_text_verbatim("PCA_TextOutput"),
                        ),
                        ui.nav_panel("M-Distance",
                            ui.output_data_frame("Ref_MDistQuartiles"),
                            ui.output_data_frame("Ref_MDistSummary"),
                            ui.input_select("ref_mdistgroup", "Average by Column:",
                                            {"Outcome":"Outcome"}, selected="Outcome"
                            ),
                            ui.output_data_frame("Ref_MDistGroup")
                        ),
                        ui.nav_panel("Normality Plot",
                            ui.input_select("ref_normalitykind", "Plot Type:", 
                                            {"pdf":"PDF", "cdf":"CDF"}, selected="pdf"
                            ),
                            ui.input_select("ref_normalitygroup", "Group by:", 
                                            {"Outcome":"Outcome"}, selected="Outcome"
                            ),
                            ui.output_plot("Ref_NormalityPlot"),
                        ),
                    ), # end ui.navset_card_pill
                ), # end ui.layout_sidebar
            ), # end ui.card
        ), # end ui.nav_panel

        ### Sample Upload
        ui.nav_panel("Sample",
            ui.card(
                ui.layout_sidebar(
                    ui.sidebar(
                    ui.input_select("samptype", "Sample File Format:", file_formats,
                        selected=default_format
                    ),  
                    ui.input_file('sampdata', 'Upload Sample Spectra',
                        multiple=False, accept=accepted_extensions
                    ),
                    ui.input_checkbox("sampbkg_check", "Apply Background?", False),
                    ui.panel_conditional("input.sampbkg_check",
                        ui.input_file('sampbkg', 'Upload Sample Background',
                            multiple=False, accept=accepted_extensions
                        )
                    ),
                    ui.input_checkbox("samp_logtrans", "Log Transform", False),
                    ui.input_checkbox("samp_trunc", "Truncation", True),
                    bg="#f8f8f8", open='always'
                    ), # end ui.sidebar
                    
                    ui.navset_card_pill(
                        ui.nav_panel("Spectra",
                            ui.output_plot("Samp_SpectraPlot"),
                            ui.output_data_frame("Samp_SpectraTable"),
                        ),
                        ui.nav_panel("M-Distance",
                            ui.output_data_frame("Samp_MDistQuartiles"),
                            ui.output_data_frame("Samp_MDistSummary"),
                            ui.input_select("samp_mdistgroup", "Average by Column:",
                                            {"Outcome":"Outcome"}, selected="Outcome"
                            ),
                            ui.output_data_frame("Samp_MDistGroup")
                        ),
                        ui.nav_panel("Normality Plot",
                            ui.input_select("samp_normalitykind", "Plot Type:",
                                            {"pdf":"PDF", "cdf":"CDF"}, selected="pdf"
                            ),
                            ui.input_select("samp_normalitygroup", "Group by:", 
                                            {"Outcome":"Outcome"}, selected="Outcome"
                            ),
                            ui.output_plot("Samp_NormalityPlot"),
                        ),
                    ), # end ui.navset_card_pill
                ), # end ui.layout_sidebar
            ), # ui.card
        ), # end ui.nav_panel
    ), # end navset_tab
) # end ui.page_fluid

### SERVER ####
def server(input: Inputs, output: Outputs, session: Session):
    ### File Upload ###
    # General Functions
    def fileRead(data_input, data_type):
        file: list[FileInfo] | None = data_input()
        if file is None:
            return None
        return read_SpectraCSV(file[0]["datapath"], format=str(data_type))

    def update_trunc_range(data_input, data_type, lower_tag, upper_tag):
        S = fileRead(data_input, data_type)
        if S is None:
            return
        ui.update_numeric(lower_tag, min=int(S.wavevals[0]), max=int(S.wavevals[-1]),
                          value=int(S.wavevals[0]))
        ui.update_numeric(upper_tag, min=int(S.wavevals[0]), max=int(S.wavevals[-1]),
                          value=int(S.wavevals[-1]))

    def get_data(data_input, data_type, bkg_check, bkg_input, logtrans_check, trunc_check):
        # Get Data and Truncate
        S = fileRead(data_input, data_type)
        if trunc_check is True and input.ref_trunc():
            S.trunc((input.ref_truncmin(), input.ref_truncmax()), True)

        # Apply Background if desired
        if bkg_check():
            bkg_file = bkg_input()
            bkg = read_SpectraCSV(bkg_file[0]["datapath"], format=str(data_type))
            bkg.trunc((input.ref_truncmin(), input.ref_truncmax()), True)
            S.apply_bkg(bkg.spc)

        # Log Transform and Output
        if logtrans_check():
            S.spc = np.log10(1/S.spc)
        return S

    # Truncation Bound Updaters
    @reactive.calc
    def Update_TruncBounds_Ref():
        S = fileRead(input.refdata, input.reftype)
        if S is None:
            return
        ui.update_numeric('ref_truncmin', min=int(S.wavevals[0]), max=int(S.wavevals[-1]),
                          value=int(S.wavevals[0]))
        ui.update_numeric('ref_truncmax', min=int(S.wavevals[0]), max=int(S.wavevals[-1]),
                          value=int(S.wavevals[-1]))
    
    # Reactive Data Get Functions
    @reactive.calc
    def Get_RefData():
        Update_TruncBounds_Ref()
        return get_data(input.refdata, input.reftype, input.refbkg_check, input.refbkg,
                        input.ref_logtrans, input.ref_trunc())

    @reactive.calc
    def Get_SampData():
        return get_data(input.sampdata, input.samptype, input.sampbkg_check, input.sampbkg,
                        input.samp_logtrans, input.samp_trunc())

    ### M-Distance ###
    @reactive.calc
    def Get_MDistModel():
        md = MDist(ncomp=input.pca_ncomp(), center=input.pca_center(),
                   scale=input.pca_scale(), ci=input.mdist_ci())
        md.fit(Get_RefData().spc)
        
        return md

    ### Plot Renders ###
    # General Functions
    def spectraplot(data_getfunc, color):
        S = data_getfunc()
        fig = plt.figure()
        plt.plot(S.wve.T, S.spc.T, linewidth=2.5, c=color)
        return fig

    def update_normalitygroups(mdsum, reforsamp='ref'):
        mdcols = dict()
        for col in mdsum.columns:
            mdcols[col] = col
        mdcols.pop("M-Distance")
        if reforsamp == 'ref':
            ui.update_select('ref_normalitygroup', label="Group by:", choices=mdcols, selected='Outcome')
        elif reforsamp == 'samp':
            ui.update_select('samp_normalitygroup', label="Group by::", choices=mdcols, selected='Outcome')
            
    def normalityplot(mdsum, plotkind, groupby, reforsamp='ref'):
        normdf = pd.DataFrame({"M-Distance":mdsum["M-Distance"], groupby:mdsum[groupby]})

        # Get boolean for cummulative parameter       
        pltkind = False
        if plotkind == 'cdf':
            pltkind = True

        # Build Facetgrid and get kdeplots
        fg = sns.FacetGrid(normdf, row=groupby, hue=groupby, aspect=15, height=0.75)
        fg.map(sns.kdeplot, 'M-Distance', bw_adjust=1, clip_on=False, fill=True, 
               alpha=0.5, linewidth=1.5, cumulative=pltkind) # kde fills
        fg.map(sns.kdeplot, 'M-Distance', bw_adjust=1, clip_on=False, color="black", 
               lw=2, cumulative=pltkind) # kde outlines
        fg.map(plt.axhline, y=0, lw=2, clip_on=False, color='black') # ridgelines

        # Labeling
        for i, ax in enumerate(fg.axes.flat):
            ax.text(-15, 0.02, normdf[groupby].to_list()[i], fontweight='bold', fontsize=15, color='black')
            ax.axvline(x=Get_MDistModel().Cutoff, linestyle='--', color='black', alpha=0.5)
        fg.set_titles("")
        fg.set(yticks=[])
        fg.despine(bottom=True, left=True)
        fg.fig.subplots_adjust(hspace=-0.3) # overlap subgraphs
        plt.setp(ax.get_xticklabels(), fontsize=15, fontweight='bold')
        plt.xlabel('M-Distance', fontweight='bold', fontsize=15)
        fg.fig.suptitle('M-Distance by ' + groupby, ha='center', fontsize=20, fontweight=20)
        
        return fg
        
    # Plot Outputs
    @render.plot
    def Ref_SpectraPlot():
        update_normalitygroups(Get_MDistModel().Summary(Get_RefData().spc, Get_RefData().var), 'ref')
        return spectraplot(Get_RefData, color='dodgerblue')
    @render.plot
    def Samp_SpectraPlot():
        update_normalitygroups(Get_MDistModel().Summary(Get_SampData().spc, Get_SampData().var), 'samp')
        return spectraplot(Get_SampData, color='darkorange')

    @render.plot
    def Ref_NormalityPlot():
        S = Get_RefData()
        mdsum = Get_MDistModel().Summary(S.spc, S.var)
        return normalityplot(mdsum, input.ref_normalitykind(), input.ref_normalitygroup(), 'ref')
    @render.plot
    def Samp_NormalityPlot():
        S = Get_SampData()
        mdsum = Get_MDistModel().Summary(S.spc, S.var)
        return normalityplot(mdsum, input.samp_normalitykind(), input.samp_normalitygroup(), 'samp')

    @render_widget 
    def PCA_ScreePlot():
        pca = Get_MDistModel().pca

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=np.arange(pca.n_components_) + 1, 
                                 y=pca.explained_variance_ratio_,
                                 line=dict(color='royalblue', width=4)
                                )
                     )
        fig.update_traces(mode='lines+markers', marker_line_width=2, marker_size=10)
        fig.update_layout(
            title=dict(text='PCA Screeplot', font=dict(size=20)), title_x=0.5,
            xaxis=dict(title=dict(text='Principal Component')),
            yaxis=dict(title=dict(text='Explained Variance (%)')),
            template='simple_white'
        )
        return fig
                      

    ### Table Renders ###
    # General Functions
    def renderDF_fromspectra(S, filters=True):
        df = S.to_df()
        df.columns = df.columns.astype('string')
        return render.DataGrid(data=df, filters=filters)

    def buildStyles(mdsum, passcol="#e6ffe6", failcol="#ffe6e6"):
        mdtest = mdsum["Outcome"].to_list()
        mdpass, mdfail = [], []
        for idx in range(len(mdtest)):
            if mdtest[idx] == 'Pass':
                mdpass.append(idx)
            elif mdtest[idx] == 'Fail':
                mdfail.append(idx)
                
        styles=[{"class": "text-center"}, # Overall Center
                {"cols": [0], "style": {"font-weight": "bold"}}, # Bold Column 1
                {"rows": mdpass, "style": {"background-color": passcol}}, # Color Passes
                {"rows": mdfail, "style": {"background-color": failcol}}  # Color Fails
               ]
        return styles
        
    def renderMDistSummary(S, filters=True, passcol="#e6ffe6", failcol="#ffe6e6", reforsamp='ref'):
        md = Get_MDistModel()
        mdsum = md.Summary(S.spc, S.var)

        # Update Groupby selection
        mdcols = dict()
        for col in mdsum.columns:
            mdcols[col] = col
        mdcols.pop("M-Distance")
        if reforsamp == 'ref':
            ui.update_select('ref_mdistgroup', label="Average by Column:", choices=mdcols, selected='Outcome')
        elif reforsamp == 'samp':
            ui.update_select('samp_mdistgroup', label="Average by Column:", choices=mdcols, selected='Outcome')
        
        styles=buildStyles(mdsum, passcol, failcol)
        return render.DataTable(data=mdsum, filters=filters, styles=styles)

    def renderMDistGroup(S, groupby='Outcome', filters=True, passcol="#e6ffe6", failcol="#ffe6e6"):
        md = Get_MDistModel()
        mdsum = md.SummaryGroupBy(S.spc, S.var, groupby)
        if groupby == "Outcome":
            mdsum.drop("Outcome", axis=1, inplace=True)
        mdsum.insert(0, groupby, mdsum.index)         
        styles=buildStyles(mdsum, passcol, failcol)
        return render.DataTable(data=mdsum, filters=filters, styles=styles)
    
    @render.data_frame
    def Ref_SpectraTable():
        return renderDF_fromspectra(Get_RefData())
    @render.data_frame  
    def Samp_SpectraTable():
        return renderDF_fromspectra(Get_SampData())

    @render.data_frame
    def Ref_MDistQuartiles():
        df = Get_MDistModel().Quartiles(Get_RefData().spc)
        df.insert(0, "Result", df.index)
        return render.DataGrid(data=df)
    @render.data_frame  
    def Samp_MDistQuartiles():
        df = Get_MDistModel().Quartiles(Get_SampData().spc)
        df.insert(0, "Result", df.index)
        return render.DataGrid(data=df)
    
    @render.data_frame
    def Ref_MDistSummary():
        return renderMDistSummary(Get_RefData(), reforsamp='ref')
    @render.data_frame  
    def Samp_MDistSummary():
        return renderMDistSummary(Get_SampData(), reforsamp='samp')

    @render.data_frame
    def Ref_MDistGroup():
        return renderMDistGroup(Get_RefData(), groupby=input.ref_mdistgroup())
    @render.data_frame  
    def Samp_MDistGroup():
        return renderMDistGroup(Get_SampData(), groupby=input.samp_mdistgroup())

app = App(app_ui, server)