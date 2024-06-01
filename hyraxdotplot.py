def index_assembly(index_file):
    """use samtools fai file to generate"""
    contig_dict = {}
    contig_list = []
    with open(index_file, 'r') as index:
        for line in index:
            row = line.rstrip("\n").split("\t")
            contig_list.append(row[0])
            contig_dict[(row[0])] = int(row[1])
    sorted_contig_list = list(dict(sorted(contig_dict.items(), key=lambda item: item[1], reverse=True)).keys())
    cumulative_length_dict = {}
    cumulative_length = 0
    for contig in sorted_contig_list:
        cumulative_length_dict[contig] = cumulative_length
        cumulative_length = cumulative_length + contig_dict[contig]
    
    #add final contig line/end line
    cumulative_length_dict["end"] = cumulative_length_dict[contig] + contig_dict[contig]
    return cumulative_length_dict


def parse_nucmer_coords_file(coords_file, threshold, x_cumulative_length_dict, y_cumulative_length_dict):    
    """parse nucmer (-T -l) coordinates file"""
    query_positions = []
    subject_positions = []
    identities = []
    strands = []
    queries = []
    subjects = []
    with open(coords_file, 'r') as f:
        for line in f:
            if line.startswith("/"): 
                pass
            elif line.startswith("NUCMER"):
                pass
            else:
                line = line.rstrip("\n")
                fields = line.split("\t")
                if len(fields) == 1:
                    pass
                elif fields[0] == "[S1]":
                    pass
                elif float(fields[6]) > float(threshold):
                    current_query_id = fields[9]
                    queries.append(current_query_id)
                    current_subject_id = fields[10]
                    subjects.append(current_subject_id)
                    query_start = int(fields[0]) + int(x_cumulative_length_dict[current_query_id])
                    query_end = int(fields[1]) + int(x_cumulative_length_dict[current_query_id])
                    subject_start = int(fields[2]) + int(y_cumulative_length_dict[current_subject_id])
                    subject_end = int(fields[3]) + int(y_cumulative_length_dict[current_subject_id])
                    identity = float(fields[6])
                    if subject_end > subject_start:    
                        strand = "+"
                    elif subject_start > subject_end:
                        strand = "-"
                    query_positions.append((query_start, query_end))
                    subject_positions.append((subject_start, subject_end))
                    identities.append(identity)
                    strands.append(strand)
    return query_positions, subject_positions, identities, strands, queries, subjects


def track(bed_file, cumulative_length_dict):
    """prepare data for loading a quantitative track, like coverage"""
    coordinates = []
    track_values = []
    contigs = []
    window_size = []
    with open(bed_file) as file:
        for line in file:
            line = line.rstrip("\n")
            fields = line.split("\t")
            start = int(fields[1]) + int(cumulative_length_dict[fields[0]])
            end = int(fields[2]) + int(cumulative_length_dict[fields[0]])
            coordinate = (start + end) // 2  # Midpoint of the window
            value = int(fields[3])
            coordinates.append(coordinate)
            track_values.append(value)
            contigs.append(fields[0])
            if len(window_size) < 1:
                window_size.append(int(end - start))
            
    track_data = {
        'position': coordinates,
        'value': track_values,
        'sequence': contigs
    }

    return track_data, window_size[0]


def plot_interactive_bokeh(coords_file, threshold, x_cumulative_length_dict, y_cumulative_length_dict, plot_title, output_file_name, plot_width, plot_height, x_feature_bed_file=None, y_feature_bed_file=None, x_track_bed_file=None, x_track_title=None, x_track_feature_name=None, x_track_colour=None, y_track_bed_file=None, y_track_title=None, y_track_feature_name=None, y_track_colour=None):
    
    """Generate interactive html dotplot"""
    
    import holoviews as hv
    from bokeh.plotting import figure, show
    from bokeh.models import ColumnDataSource, Span, HoverTool, NumeralTickFormatter, TapTool, CustomJS, Button, BoxAnnotation, Div, Select
    from bokeh.layouts import column, row
    from matplotlib.colors import Normalize, LinearSegmentedColormap
    import matplotlib.pyplot as plt
    from bokeh.io import output_file
    
    query_positions, subject_positions, identities, strands, queries, subjects = parse_nucmer_coords_file(coords_file, threshold, x_cumulative_length_dict, y_cumulative_length_dict)
    hv.extension('bokeh')
    p = figure(title="Interactive Dot Plot", width=plot_width, height=plot_height, min_border_left=100, x_axis_label='Position (bp)', y_axis_label='Position (bp)', tools="pan, reset, box_zoom")
    #define color palette
    colors = ['orange', 'red', 'black']
    markers = [0, 0.5, 1]
    norm = plt.Normalize(vmin=min(identities), vmax=max(identities))
    blended_cmap = LinearSegmentedColormap.from_list("blended_cmap", list(zip(markers, colors)))
    
    #remove bokeh grids
    p.xgrid.visible = False
    p.ygrid.visible = False
    for x_position in x_cumulative_length_dict.values():
        vline = Span(location=x_position, dimension='height', line_color='grey', line_width=0.5)
        p.add_layout(vline)       
    for y_position in y_cumulative_length_dict.values():
        hline = Span(location=y_position, dimension='width', line_color='grey', line_width=0.5)
        p.add_layout(hline) 
    
    #plot dot plot as multi-line bokeh plot
    x = []
    y = []
    colors = []
    query_contigs = []
    subject_contigs = []
    for (q_start, q_end), (s_start, s_end), identity, strand, query, subject in zip(query_positions, subject_positions, identities, strands, queries, subjects):
        x.append([q_start, q_end])
        if strand == "+":
            y.append([s_start, s_end])
        elif strand == "-":
            y.append([s_start, s_start + (q_start - q_end)])
        

        color=blended_cmap(norm(identity))
        hex_color = "#{:02X}{:02X}{:02X}".format(int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))
        colors.append(hex_color)

        query_contigs.append(query)
        subject_contigs.append(subject)
    
    source = ColumnDataSource(data={'x': x, 'y': y, 'colors': colors, 'query': query_contigs, 'subject': subject_contigs, 'identity':identities})
    multi_line = p.multi_line('x', 'y', color='colors', source=source, line_width=1)

    # Add the HoverTool to the plot
    hover_lines = HoverTool(renderers=[multi_line], tooltips=[("x", "@query"), ("y", "@subject"), ("identity", "@identity")])
    p.add_tools(hover_lines)

    # Configure the x-axis and y-axis tick formatters
    p.xaxis.formatter = NumeralTickFormatter(format="0")
    p.yaxis.formatter = NumeralTickFormatter(format="0")
    p.yaxis.formatter = NumeralTickFormatter(format="0,0")
    p.xaxis.formatter = NumeralTickFormatter(format="0,0")

    # Plot repeats
    def plot_feature(bed_file, cumulative_length_dict, axis, color):
        features = []
        with open(bed_file) as file:
            for line in file:
                line = line.rstrip("\n")
                fields = line.split("\t")
                feature_start = int(fields[1]) + int(cumulative_length_dict[fields[0]])
                feature_end = int(fields[2]) + int(cumulative_length_dict[fields[0]])
                if axis == "x":
                    box = BoxAnnotation(left=feature_start, right=feature_end, fill_color=color, fill_alpha=0.3)
                elif axis == "y":
                    box = BoxAnnotation(bottom=feature_start, top=feature_end, fill_color=color, fill_alpha=0.3)
                features.append(box)
                p.add_layout(box)
        return features

    # Initial feature color
    initial_color = 'green'
    x_features = []
    y_features = []
    if x_feature_bed_file:
        x_features = plot_feature(x_feature_bed_file, x_cumulative_length_dict, "x", initial_color)
    if y_feature_bed_file:
        y_features = plot_feature(y_feature_bed_file, y_cumulative_length_dict, "y", initial_color)


    # Div widget to display selected contigs
    selected_contigs_div = Div(text="<b>Selected Contigs:</b><br>", width=200, height=400)

    # JavaScript callback to update Div with selected contigs
    callback = CustomJS(args=dict(source=source, contigs_div=selected_contigs_div), code="""
        var indices = source.selected.indices;
        var data = source.data;
        var selected_contigs = [];
        
        // Retrieve previously selected contigs
        var existing_contigs = contigs_div.text.split('<br>').slice(1);
        for (var i = 0; i < existing_contigs.length; i++) {
            if (existing_contigs[i] !== "") {
                selected_contigs.push(existing_contigs[i]);
            }
        }

        // Toggle selection of newly clicked contigs
        for (var i = 0; i < indices.length; i++) {
            var contig = data['query'][indices[i]];
            var contigIndex = selected_contigs.indexOf(contig);
            if (contigIndex === -1) {
                selected_contigs.push(contig); // Add contig if not already selected
            } else {
                selected_contigs.splice(contigIndex, 1); // Remove contig if already selected
            }
        }

        // Update the Div with selected contigs
        contigs_div.text = "<b>Selected Contigs:</b><br>" + selected_contigs.join('<br>');
    """)

    # Add TapTool for selection with the callback
    tap_tool = TapTool(callback=callback)
    p.add_tools(tap_tool)

    # JavaScript callback to download selected contigs
    download_callback = CustomJS(args=dict(contigs_div=selected_contigs_div), code="""
        // Get the list of selected contigs from the Div
        var selected_contigs = contigs_div.text.split('<br>').slice(1);
        if (selected_contigs.length === 0) {
            alert("No contigs selected");
            return;
        }

        // Create a downloadable file
        var filetext = selected_contigs.join('\\n');
        var filename = 'selected_contigs.txt';

        var blob = new Blob([filetext], { type: 'text/plain;charset=utf-8;' });
        if (navigator.msSaveBlob) { // IE 10+
            navigator.msSaveBlob(blob, filename);
        } else {
            var link = document.createElement("a");
            if (link.download !== undefined) {
                var url = URL.createObjectURL(blob);
                link.setAttribute("href", url);
                link.setAttribute("download", filename);
                link.style.visibility = 'hidden';
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
            }
        }
    """)

    # Add a button to download selected contigs
    download_button = Button(label="Download Selected Contigs", button_type="success")
    download_button.js_on_click(download_callback)


    if x_feature_bed_file:
        # JavaScript callback to toggle feature visibility
        x_toggle_callback = CustomJS(args=dict(features=x_features), code="""
            for (var i = 0; i < features.length; i++) {
                features[i].visible = !features[i].visible;
            }
        """)

        # Add a button to toggle feature visibility
        x_toggle_button = Button(label="Toggle X Bed File", button_type="warning")
        x_toggle_button.js_on_click(x_toggle_callback)

        # JavaScript callback to update feature colors
        x_update_color_callback = CustomJS(args=dict(features=x_features), code="""
            var color = cb_obj.value;
            for (var i = 0; i < features.length; i++) {
                features[i].fill_color = color;
            }
        """)

        # Add a dropdown menu to select feature color
        x_color_select = Select(title="X Bed File Colour:", value=initial_color, options=["green", "blue", "red", "purple", "yellow"])
        x_color_select.js_on_change("value", x_update_color_callback)
    
    if y_feature_bed_file:
        # JavaScript callback to toggle feature visibility
        y_toggle_callback = CustomJS(args=dict(features=y_features), code="""
            for (var i = 0; i < features.length; i++) {
                features[i].visible = !features[i].visible;
            }
        """)

        # Add a button to toggle feature visibility
        y_toggle_button = Button(label="Toggle Y Bed File", button_type="warning")
        y_toggle_button.js_on_click(y_toggle_callback)

        # JavaScript callback to update feature colors
        y_update_color_callback = CustomJS(args=dict(features=y_features), code="""
            var color = cb_obj.value;
            for (var i = 0; i < features.length; i++) {
                features[i].fill_color = color;
            }
        """)

        # Add a dropdown menu to select feature color
        y_color_select = Select(title="Y Bed File Colour:", value=initial_color, options=["green", "blue", "red", "purple", "yellow"])
        y_color_select.js_on_change("value", y_update_color_callback)
    
    if x_track_bed_file:
        # Generate x track data
        x_track_data, x_window_size = track(x_track_bed_file, x_cumulative_length_dict)
        source = ColumnDataSource(data=x_track_data)

        # Create a Bokeh plot
        j = figure(title=x_track_title, y_axis_label=x_track_feature_name, width=plot_width, height=150, x_range=p.x_range, min_border_left=100, tools="pan")
        # Add vertical bars for the quantitative track
        vbar = j.vbar(x='position', top='value', width=x_window_size, source=source, line_color=x_track_colour, fill_color=x_track_colour)
        # Add the HoverTool to the plot
        hover = HoverTool(renderers=[vbar], tooltips=[("Sequence", "@sequence"), ("Position", "@position"), (x_track_feature_name, "@value")])
        j.add_tools(hover)

        # Configure the x-axis and y-axis tick formatters
        j.xaxis.formatter = NumeralTickFormatter(format="0")
        j.yaxis.formatter = NumeralTickFormatter(format="0")
        j.yaxis.formatter = NumeralTickFormatter(format="0,0") #add commas to numbers

        #remove shared axis ticks
        j.xaxis.major_tick_line_color = None
        j.xaxis.major_label_text_font_size = '0pt'

        #remove bokeh grids and add grid lines to match contig boundaries
        j.xgrid.visible = False
        j.ygrid.visible = True
        j.ygrid.grid_line_dash = [4, 4]
        for x_position in x_cumulative_length_dict.values():
            vline = Span(location=x_position, dimension='height', line_color='grey', line_width=0.5)
            j.add_layout(vline) 
    
    if y_track_bed_file:
        # Generate y track data
        y_track_data, y_window_size = track(y_track_bed_file, y_cumulative_length_dict)
        source = ColumnDataSource(data=y_track_data)

        # Create a Bokeh plot
        c = figure(title=y_track_title, x_axis_label=y_track_feature_name, height=plot_height, width=150, y_range=p.y_range, tools="pan")
        # Add vertical bars for the y track
        hbar = c.hbar(y='position', right='value', height=y_window_size, source=source, line_color=y_track_colour, fill_color=y_track_colour)
        # Add the HoverTool to the plot
        hover = HoverTool(renderers=[hbar], tooltips=[("Sequence", "@sequence"), ("Position", "@position"), (y_track_feature_name, "@value")])
        c.add_tools(hover)

        # Configure the x-axis and y-axis tick formatters
        c.xaxis.formatter = NumeralTickFormatter(format="0")
        c.yaxis.formatter = NumeralTickFormatter(format="0")
        c.xaxis.formatter = NumeralTickFormatter(format="0,0")

        #remove shared axis ticks
        c.yaxis.major_tick_line_color = None
        c.yaxis.major_label_text_font_size = '0pt'
        c.xaxis.major_label_text_font_size = '6pt'

        #remove bokeh grids
        c.xgrid.visible = True
        c.ygrid.visible = False
        c.xgrid.grid_line_dash = [4, 4]
        for y_position in y_cumulative_length_dict.values():
            hline = Span(location=y_position, dimension='width', line_color='grey', line_width=0.5)
            c.add_layout(hline) 



    title = Div(text=f"<h1>{plot_title}</h1>", width=plot_width, height=50, styles={'text-align': 'center'})
    
    #generate final layout
    if x_feature_bed_file and y_feature_bed_file:
        if y_track_bed_file and x_track_bed_file:
            layout = column(title,
                        j,
                        row(p, c), 
                        row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
        elif y_track_bed_file:
            layout = column(title,
                        row(p, c), 
                        row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
        elif x_track_bed_file:
            layout = column(title,
                        j,
                        p,
                        row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
        else:
            layout = column(title,
                        p,
                        row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
    elif x_feature_bed_file:
        if y_track_bed_file and x_track_bed_file:
            layout = column(title,
                        j,
                        row(p, c), 
                        row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
        elif y_track_bed_file:
            layout = column(title,
                        row(p, c), 
                        row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
        elif x_track_bed_file:
            layout = column(title,
                        j,
                        p,
                        row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
        else:
            layout = column(title,
                        p,
                        row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
            
    elif y_feature_bed_file:
        if y_track_bed_file and x_track_bed_file:
            layout = column(title,
                        j,
                        row(p, c), 
                        row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
        elif y_track_bed_file:
            layout = column(title,
                        row(p, c), 
                        row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
        elif x_track_bed_file:
            layout = column(title,
                        j,
                        p,
                        row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
        else:
            layout = column(title,
                        p,
                        row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
    
    else:
        if y_track_bed_file and x_track_bed_file:
            layout = column(title,
                        j,
                        row(p, c), 
                        row(selected_contigs_div, download_button))
        elif y_track_bed_file:
            layout = column(title,
                        row(p, c), 
                        row(selected_contigs_div, download_button))
        elif x_track_bed_file:
            layout = column(title,
                        j,
                        p,
                        row(selected_contigs_div, download_button))
        else:
            layout = column(title,
                        p,
                        row(selected_contigs_div, download_button))


    # Save and show final plot
    output_file(output_file_name)
    
    #show(layout)


def main():
    
    import argparse
    import sys


    print("""
    █░█ █▄█ █▀█ ▄▀█ ▀▄▀ █▀▄ █▀█ ▀█▀ █▀█ █░░ █▀█ ▀█▀
    █▀█ ░█░ █▀▄ █▀█ █░█ █▄▀ █▄█ ░█░ █▀▀ █▄▄ █▄█ ░█░
    """)
    
    
    print("""
    easily generate interactive dot plots, with em-
    bedded feature bed files that can be loaded as 
    highlighted regions; and quantitative bed files 
    which can be loaded as tracks (e.g. read depth).
          
    Please contact ak37@sanger.ac.uk for issues, 
    questions and requests.
    ________________________________________________
                
                  version 1.0.0
    """)

    print("""
             ▒▒▒▒▒▒░▒▒▒▒▒▒                 
            ▒▒▒▒▒▒▒▒▒▒▒▒▒░▒▒▒▒▒▒▒▒▒▓▓     
          ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒     
          ▓▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▓▓▒▒▒▒▒▒   
          ▓▒▒▒▓▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▓▓▓▓▒▒▒▒▓▓▓ 
          ▓▓▓▓▓▓▓▒▓▒▒▒▒▒▒▒▒▓▓▓▓▒▒▒▒▓▓██▓▓▓▓
          ▓▒▒▓▓▓▓▒▓▓▓▓▓▓▓▓▒▒▓▓▓▓▒▓▒▒▓▓▓▓▓▓█
           ▓▓▓▓▓▓▓▒▓▓▓▓▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓███
             ▓▓▓██▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓          
             ▓▓▓▓   ▓▓▓▓ ▓▓█▓█ ▓▓          
              ▓▓▓▓▓        ▓▓▓ ▓▓▓▓        
                ▓▓         ▓▒▓▓            
    """)

    parser = argparse.ArgumentParser()
    parser.add_argument("--coords_file", metavar="FILE", help="nucmer coordinates file generated with -T and -l flags; x-axis fasta should be the query sequence and y-axis fasta should be the subject sequence (i.e. nucmer x.fasta y.fasta)")
    parser.add_argument("--x_index_file", metavar="FILE", help="samtools .fai file of fasta to be on x-axis of dotplot")
    parser.add_argument("--y_index_file", metavar="FILE", help="samtools .fai file of fasta to be on y-axis of dotplot")
    parser.add_argument("--threshold", metavar="FLOAT", type=float, default=90, help="(optional) minimum nucleotide identity of matches to be plotted; default is 90")
    parser.add_argument("--plot_title", metavar="STR", default="HyraxDotPlot", help="(optional) plot title; default is HyraxDotPlot")
    parser.add_argument("--output", metavar="STR", default="hyraxdotplot.html", help="(optional) output file name; default is hyraxdotplot.html")
    parser.add_argument("--x_feature_bed_file", metavar="FILE", help="(optional) qualitative bed file for x-axis fasta to be loaded as highlighted regions on dot plot")
    parser.add_argument("--y_feature_bed_file", metavar="FILE", help="(optional) qualitative bed file for y-axis fasta to be loaded as highlighted regions on dot plot")
    parser.add_argument("--plot_width", metavar="INT", type=int, default=800, help="(optional) plot width; default is 800")
    parser.add_argument("--plot_height", metavar="INT", type=int, default=600, help="(optional) plot height; default is 600")
    parser.add_argument("--x_track_bed_file", metavar="FILE", help="quantitative bed file to be loaded as a track for x-axis fasta, e.g. windowed coverage bed file")
    parser.add_argument("--x_track_title", metavar="STR", default="Feature Plot", help="(optional) plot title for track loaded on x-axis fasta")
    parser.add_argument("--x_track_feature_name", metavar="STR", default="Feature", help="(optional) feature name for track loaded on x-axis fasta (i.e. what hover tool will call feature on track)")
    parser.add_argument("--x_track_colour", metavar="STR", default="#56B4E9", help="(optional) colour of track loaded on x-axis fasta; default is #56B4E9")
    parser.add_argument("--y_track_bed_file", metavar="FILE", help="quantitative bed file to be loaded as a track for y-axis fasta, e.g. windowed coverage bed file")
    parser.add_argument("--y_track_title", metavar="STR", default="Feature Plot", help="(optional) plot title for track loaded on y-axis fasta")
    parser.add_argument("--y_track_feature_name", metavar="STR", default="Feature", help="(optional) feature name for track loaded on y-axis fasta (i.e. what hover tool will call feature on track)")
    parser.add_argument("--y_track_colour", metavar="STR", default="#f49ac2", help="(optional) colour of track loaded on y-axis fasta; default is #f49ac2")

    args = parser.parse_args()

    # Check if positional arguments are missing and print an error message if they are
    missing_args = []
    if not args.coords_file:
        missing_args.append("--coords_file")
    if not args.x_index_file:
        missing_args.append("--x_index_file")
    if not args.y_index_file:
        missing_args.append("--y_index_file")

    if missing_args:
        print(f"Error: Missing required positional arguments: {', '.join(missing_args)}")
        parser.print_help()
        sys.exit(1)
    
    xcumulative_length_dict = index_assembly(args.x_index_file)
    ycumulative_length_dict = index_assembly(args.y_index_file)

    plot_interactive_bokeh(coords_file=args.coords_file, threshold=args.threshold, x_cumulative_length_dict=xcumulative_length_dict, y_cumulative_length_dict=ycumulative_length_dict, plot_title=args.plot_title, output_file_name=args.output, plot_width=args.plot_width, plot_height=args.plot_height, x_feature_bed_file=args.x_feature_bed_file, y_feature_bed_file=args.y_feature_bed_file, x_track_bed_file=args.x_track_bed_file, x_track_title=args.x_track_title, x_track_feature_name=args.x_track_feature_name, x_track_colour=args.x_track_colour, y_track_bed_file=args.y_track_bed_file, y_track_title=args.y_track_title, y_track_feature_name=args.y_track_feature_name, y_track_colour=args.y_track_colour)

    
    print("""Plotting finished, enjoy your plot""")

if __name__ == "__main__":
    main()
