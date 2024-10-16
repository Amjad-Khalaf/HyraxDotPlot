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

def parse_nucmer_coords_file(coords_file, threshold, size_threshold, x_cumulative_length_dict, y_cumulative_length_dict):    
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
                elif float(fields[6]) > float(threshold) and fields[9] in x_cumulative_length_dict.keys() and fields[10] in y_cumulative_length_dict.keys() and abs(int(fields[0]) - int(fields[1])) >= size_threshold:
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

def parse_paf_file(paf_file, threshold, size_threshold, x_cumulative_length_dict, y_cumulative_length_dict):    
    """parse paf file"""
    query_positions = []
    subject_positions = []
    identities = []
    strands = []
    queries = []
    subjects = []
    with open(paf_file, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            fields = line.split("\t")
            for i in range(0, len(fields)):
                if "de" in fields[i]:
                    identity = 100 - (float(fields[i].split(":")[2])*100)
            if identity > float(threshold) and fields[0] in x_cumulative_length_dict.keys() and fields[5] in y_cumulative_length_dict.keys() and abs(int(fields[2]) - int(fields[3])) >= size_threshold:
                current_query_id = fields[0]
                queries.append(current_query_id)
                current_subject_id = fields[5]
                subjects.append(current_subject_id)
                query_start = int(fields[2]) + int(x_cumulative_length_dict[current_query_id])
                query_end = int(fields[3]) + int(x_cumulative_length_dict[current_query_id])
                subject_start = int(fields[7]) + int(y_cumulative_length_dict[current_subject_id])
                subject_end = int(fields[8]) + int(y_cumulative_length_dict[current_subject_id])
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
            if fields[0] in cumulative_length_dict.keys():
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

def annotation(bed_file, cumulative_length_dict):
    """prepare data for loading an annotation track"""
    start_positions = []
    end_positions = []
    contigs = []
    names = []
    strands = []
    object_centres = []
    object_y_coordinates = []
    colors = []
    with open(bed_file) as file:
        for line in file:
            line = line.rstrip("\n")
            fields = line.split("\t")
            if fields[0] in cumulative_length_dict.keys():
                start = int(fields[1]) + int(cumulative_length_dict[fields[0]])
                end = int(fields[2]) + int(cumulative_length_dict[fields[0]])
                start_positions.append(start)
                end_positions.append(end)
                name = fields[3]
                names.append(name)
                strand = fields[5]
                strands.append(strand)
                contigs.append(fields[0])
                object_centres.append([start, end])
                object_y_coordinates.append([0,0])
                if strand == '+':
                    colors.append('#EFA647')
                elif strand == '-':
                    colors.append('#5A70A3')
    # Create ColumnDataSource with all data
    x_annotations_data={
        'x': object_centres,
        'y': object_y_coordinates,
        'color': colors,
        'name': names,
        'start': start_positions,
        'end': end_positions,
        'strand': strands,
        'location': contigs
        }      
    return x_annotations_data, start_positions, end_positions, strands

def plot_interactive_bokeh(threshold, size_threshold, x_cumulative_length_dict, y_cumulative_length_dict, plot_title, output_prefix, plot_width, plot_height, coords_file=None, paf_file=None, x_feature_bed_file=None, y_feature_bed_file=None, x_track_bed_file=None, x_track_title=None, x_track_feature_name=None, x_track_colour=None, y_track_bed_file=None, y_track_title=None, y_track_feature_name=None, y_track_colour=None, x_annotation_bed_file=None, y_annotation_bed_file=None, curation_mode=False):
    
    """Generate interactive html dotplot"""
    
    import holoviews as hv
    from bokeh.plotting import figure, show, save
    from bokeh.models import HoverTool, FixedTicker, Label, ColumnDataSource, Span, NumeralTickFormatter, TapTool, CustomJS, Button, BoxAnnotation, Div, Select, Range1d, WheelZoomTool, CrosshairTool, LinearColorMapper, ColorBar
    from bokeh.layouts import column, row
    from matplotlib.colors import Normalize, LinearSegmentedColormap
    import matplotlib.pyplot as plt
    from bokeh.transform import linear_cmap
    from bokeh.io import output_file, export_png, export_svg
    
    if coords_file:
        query_positions, subject_positions, identities, strands, queries, subjects = parse_nucmer_coords_file(coords_file, threshold, size_threshold, x_cumulative_length_dict, y_cumulative_length_dict)
    if paf_file:
        query_positions, subject_positions, identities, strands, queries, subjects = parse_paf_file(paf_file, threshold, size_threshold, x_cumulative_length_dict, y_cumulative_length_dict)
    
    hv.extension('bokeh')
    
    # Create the Bokeh figure
    p = figure(width=plot_width, height=plot_height, min_border_left=100, x_axis_label='Position (bp)', y_axis_label='Position (bp)', tools="pan, reset, box_zoom, wheel_zoom, crosshair")
    p.toolbar.active_scroll = p.select_one(WheelZoomTool)  # Enable wheel zoom as the active scroll tool
    p.output_backend="svg"
    # Define color palette and normalization
    colors = ['orange', 'red', 'black']
    markers = [0, 0.5, 1]
    norm = plt.Normalize(vmin=min(identities), vmax=max(identities))
    blended_cmap = LinearSegmentedColormap.from_list("blended_cmap", list(zip(markers, colors)))

    # Convert the Matplotlib colormap to a list of hex colors
    cmap_hex = [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}" for r, g, b, a in blended_cmap(range(blended_cmap.N))]

    # Remove grid lines
    p.xgrid.visible = False
    p.ygrid.visible = False

    # Add vertical and horizontal lines for contig boundaries
    for x_position in x_cumulative_length_dict.values():
        vline = Span(location=x_position, dimension='height', line_color='grey', line_width=0.5)
        p.add_layout(vline)
    for y_position in y_cumulative_length_dict.values():
        hline = Span(location=y_position, dimension='width', line_color='grey', line_width=0.5)
        p.add_layout(hline)

    # Prepare data for plotting
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
        
        # Map identity values to colors using the custom color map
        color = blended_cmap(norm(identity))
        hex_color = "#{:02X}{:02X}{:02X}".format(int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))
        colors.append(hex_color)

        query_contigs.append(query)
        subject_contigs.append(subject)

    # Create a LinearColorMapper for the color bar and line coloring
    mapper = LinearColorMapper(palette=cmap_hex, low=min(identities), high=max(identities))

    # Create a ColumnDataSource for the plot
    source = ColumnDataSource(data={'x': x, 'y': y, 'colors': colors, 'query': query_contigs, 'subject': subject_contigs, 'identity': identities})
    
    # Add multi-line plot to the figure
    multi_line = p.multi_line('x', 'y', color='colors', source=source, line_width=1)

    # Add the HoverTool to the plot
    hover_lines = HoverTool(renderers=[multi_line], tooltips=[("x", "@query"), ("y", "@subject"), ("identity", "@identity")])
    p.add_tools(hover_lines)

    # Configure the x-axis and y-axis tick formatters
    p.xaxis.formatter = NumeralTickFormatter(format="0,0")
    p.yaxis.formatter = NumeralTickFormatter(format="0,0")

    # Add autohide toolbar
    p.toolbar.autohide = True

    # Add a color bar to the plot
    color_bar_width = int(plot_width*0.8)
    color_bar = ColorBar(color_mapper=mapper, width=color_bar_width, location=(0, 0), title="Nucleotide Identity (%)", bar_line_color='black', major_tick_line_color='black')
    p.add_layout(color_bar, 'below')

    # Plot features
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

    if curation_mode == True:
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
        y_axis_maximum_value = max(x_track_data["value"])
        j = figure(title=x_track_title, y_axis_label=x_track_feature_name, width=plot_width, height=150, x_range=p.x_range, min_border_left=100, tools="pan", y_axis_type="log",  y_range=[1,y_axis_maximum_value])
        j.output_backend="svg"
        # Add vertical bars for the quantitative track
        vbar = j.vbar(x='position', top='value', width=x_window_size, source=source, line_color=x_track_colour, fill_color=x_track_colour, bottom=1e-10)
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

        #add autohide
        j.toolbar.autohide = True
    
    if y_track_bed_file:
        # Generate y track data
        y_track_data, y_window_size = track(y_track_bed_file, y_cumulative_length_dict)
        source = ColumnDataSource(data=y_track_data)

        # Create a Bokeh plot
        x_axis_maximum_value = max(y_track_data["value"])
        c = figure(title=y_track_title, x_axis_label=y_track_feature_name, height=plot_height, width=200, y_range=p.y_range, tools="pan", x_axis_type="log", x_range=[1, x_axis_maximum_value])
        c.output_backend="svg"
        # Add vertical bars for the y track
        hbar = c.hbar(y='position', right='value', height=y_window_size, source=source, line_color=y_track_colour, fill_color=y_track_colour, left=1e-10)
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
        
        #add autohide
        c.toolbar.autohide = True

    if x_annotation_bed_file:
        # Generate x annotation data
        x_annotations_data, start_positions, end_positions, strands= annotation(x_annotation_bed_file, x_cumulative_length_dict)
        source = ColumnDataSource(data=x_annotations_data)

        a_x = figure(title="Annotation Track", width=plot_width, height=100, x_range=p.x_range, min_border_left=100, tools="pan", y_range=Range1d(-0.5, 0.5))
        multi_line = a_x.multi_line('x', 'y', color='color', source=source, line_width=25)
        a_x.output_backend="svg"

        # Add the HoverTool to the plot
        hover_lines = HoverTool(renderers=[multi_line], tooltips=[("Name", "@name"), ("Start Position", "@start"), ("End Position", "@end"), ("Strand", "@strand"), ("Location", "@location")])
        a_x.add_tools(hover_lines)

        #remove axes ticks
        a_x.xaxis.major_tick_line_color = None
        a_x.xaxis.major_label_text_font_size = '0pt'
        a_x.yaxis.major_tick_line_color = None
        a_x.yaxis.major_label_text_font_size = '0pt'
        
        #remove bokeh grids and add grid lines to match contig boundaries
        a_x.xgrid.visible = False
        a_x.ygrid.visible = True
        a_x.ygrid.grid_line_dash = [4, 4]
        for x_position in x_cumulative_length_dict.values():
            vline = Span(location=x_position, dimension='height', line_color='grey', line_width=0.5)
            a_x.add_layout(vline) 
        
        #add autohide
        a_x.toolbar.autohide = True

    if y_annotation_bed_file:
        # Generate x annotation data
        y_annotations_data, start_positions, end_positions, strands= annotation(y_annotation_bed_file, y_cumulative_length_dict)
        source = ColumnDataSource(data=y_annotations_data)

        a_y = figure(title="Annotation Track", height=plot_height, width=150, y_range=p.y_range, tools="pan", x_range=Range1d(-0.5, 0.5))
        multi_line = a_y.multi_line('y', 'x', color='color', source=source, line_width=25)
        a_y.output_backend="svg"

        # Add the HoverTool to the plot
        hover_lines = HoverTool(renderers=[multi_line], tooltips=[("Name", "@name"), ("Start Position", "@start"), ("End Position", "@end"), ("Strand", "@strand"), ("Location", "@location")])
        a_y.add_tools(hover_lines)

        #remove axes ticks
        a_y.xaxis.major_tick_line_color = None
        a_y.xaxis.major_label_text_font_size = '0pt'
        a_y.yaxis.major_tick_line_color = None
        a_y.yaxis.major_label_text_font_size = '0pt'
        
        #remove bokeh grids and add grid lines to match contig boundaries
        a_y.xgrid.visible = True
        a_y.ygrid.visible = False
        a_y.xgrid.grid_line_dash = [4, 4]
        for y_position in y_cumulative_length_dict.values():
            hline = Span(location=y_position, dimension='width', line_color='grey', line_width=0.5)
            a_y.add_layout(hline) 
    
        #add autohide
        a_y.toolbar.autohide = True

    title = Div(text=f"<h1>{plot_title}</h1>", width=plot_width, height=50, styles={'text-align': 'center'})
    
    #generate final layout
    
    if curation_mode == True:
        if x_feature_bed_file and y_feature_bed_file:
            
            if y_track_bed_file and x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, j, a_x, row(p, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, a_y, c))
                    else:
                        layout = column(title, j, row(p, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, c))
                
            elif y_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, a_x, row(p, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, row(p, a_y, c))
                    else:
                        layout = column(title, row(p, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, row(p, c))
            
            elif x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, a_y))
                    else:
                        layout = column(title, j, a_x, p, row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, a_y))
                    else:
                        layout = column(title, j, p, row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, p)
            
            else:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, a_y))
                    else:
                        layout = column(title, a_x, p, row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, row(p, a_y))
                    else:
                        layout = column(title, p, row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, p)
        
        elif x_feature_bed_file:
            if y_track_bed_file and x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y, c), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, j, a_x, row(p, c), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y, c), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, a_y, c))
                    else:
                        layout = column(title, j, row(p, c), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, c))

            elif y_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y, c), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, a_x, row(p, c), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y, c), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, row(p, a_y, c))
                    else:
                        layout = column(title, row(p, c), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, row(p, c))
                
            elif x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, a_y))
                    else:
                        layout = column(title, j, a_x, p, row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, a_y))
                    else:
                        layout = column(title, j, p, row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, p)
            
            else:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, a_y))
                    else:
                        layout = column(title, a_x, p, row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y), row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, row(p, a_y))
                    else:
                        layout = column(title, p, row(column(x_toggle_button, x_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, p)
                
        elif y_feature_bed_file:
            if y_track_bed_file and x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y, c), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, j, a_x, row(p, c), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y, c), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, a_y, c))
                    else:
                        layout = column(title, j, row(p, c), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, c))
            
            elif y_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y, c), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, a_x, row(p, c), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y, c), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, row(p, a_y, c))
                    else:
                        layout = column(title, row(p, c), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, row(p, c))
            
            elif x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, a_y))
                    else:
                        layout = column(title, j, a_x, p, row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, a_y))
                    else:
                        layout = column(title, j, p, row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, j, p)
            
            else:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, a_y))
                    else:
                        layout = column(title, a_x, p, row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y), row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, row(p, a_y))
                    else:
                        layout = column(title, p, row(column(y_toggle_button, y_color_select), selected_contigs_div, download_button))
                        static_layout = column(title, p)
        
        else:
            if y_track_bed_file and x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y, c), row(selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, j, a_x, row(p, c), row(selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y, c), row(selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, a_y, c))
                    else:
                        layout = column(title, j, row(p, c), row(selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, c))
            
            elif y_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y, c), row(selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, a_x, row(p, c), row(selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y, c), row(selected_contigs_div, download_button))
                        static_layout = column(title, row(p, a_y, c))
                    else:
                        layout = column(title, row(p, c), row(selected_contigs_div, download_button))
                        static_layout = column(title, row(p, c))

            elif x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y), row(selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, row(p, a_y))
                    else:
                        layout = column(title, j, a_x, p, row(selected_contigs_div, download_button))
                        static_layout = column(title, j, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y), row(selected_contigs_div, download_button))
                        static_layout = column(title, j, row(p, a_y))
                    else:
                        layout = column(title, j, p, row(selected_contigs_div, download_button))
                        static_layout = column(title, j, p)
            
            else:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y), row(selected_contigs_div, download_button))
                        static_layout = column(title, a_x, row(p, a_y))
                    else:
                        layout = column(title, a_x, p, row(selected_contigs_div, download_button))
                        static_layout = column(title, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y), row(selected_contigs_div, download_button))
                        static_layout = column(title, row(p, a_y))
                    else:
                        layout = column(title, p, row(selected_contigs_div, download_button))
                        static_layout = column(title, p)
    
    elif curation_mode == False:
        if x_feature_bed_file and y_feature_bed_file:
            
            if y_track_bed_file and x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, j, a_x, row(p, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, row(p, a_y, c))
                    else:
                        layout = column(title, j, row(p, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, row(p, c))
                
            elif y_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, a_x, row(p, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, row(p, a_y, c))
                    else:
                        layout = column(title, row(p, c), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, row(p, c))
            
            elif x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, a_x, row(p, a_y))
                    else:
                        layout = column(title, j, a_x, p, row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, row(p, a_y))
                    else:
                        layout = column(title, j, p, row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, p)
            
            else:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, a_x, row(p, a_y))
                    else:
                        layout = column(title, a_x, p, row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y), row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, row(p, a_y))
                    else:
                        layout = column(title, p, row(column(x_toggle_button, x_color_select), column(y_toggle_button, y_color_select)))
                        static_layout = column(title, p)
        
        elif x_feature_bed_file:
            if y_track_bed_file and x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y, c), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, j, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, j, a_x, row(p, c), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, j, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y, c), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, j, row(p, a_y, c))
                    else:
                        layout = column(title, j, row(p, c), row(column(x_toggle_button, x_color_select)))

            elif y_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y, c), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, a_x, row(p, c), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y, c), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, row(p, a_y, c))
                    else:
                        layout = column(title, row(p, c), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, row(p, c))
                
            elif x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, j, a_x, row(p, a_y))
                    else:
                        layout = column(title, j, a_x, p, row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, j, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, j, row(p, a_y))
                    else:
                        layout = column(title, j, p, row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, j, p)
            
            else:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, a_x, row(p, a_y))
                    else:
                        layout = column(title, a_x, p, row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y), row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, row(p, a_y))
                    else:
                        layout = column(title, p, row(column(x_toggle_button, x_color_select)))
                        static_layout = column(title, p)
                
        elif y_feature_bed_file:
            if y_track_bed_file and x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y, c), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, j, a_x, row(p, c), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y, c), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, row(p, a_y, c))
                    else:
                        layout = column(title, j, row(p, c), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, row(p, c))
            
            elif y_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y, c), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, a_x, row(p, c), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y, c), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, row(p, a_y, c))
                    else:
                        layout = column(title, row(p, c), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, row(p, c))
            
            elif x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, a_x, row(p, a_y))
                    else:
                        layout = column(title, j, a_x, p, row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, row(p, a_y))
                    else:
                        layout = column(title, j, p, row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, j, p)
            
            else:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, a_x, row(p, a_y))
                    else:
                        layout = column(title, a_x, p, row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y), row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, row(p, a_y))
                    else:
                        layout = column(title, p, row(column(y_toggle_button, y_color_select)))
                        static_layout = column(title, p)
        
        else:
            if y_track_bed_file and x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y, c))
                        static_layout = column(title, j, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, j, a_x, row(p, c))
                        static_layout = column(title, j, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y, c))
                        static_layout = column(title, j, row(p, a_y, c))
                    else:
                        layout = column(title, j, row(p, c))
                        static_layout = column(title, j, row(p, c))
            
            elif y_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y, c))
                        static_layout = column(title, a_x, row(p, a_y, c))
                    else:
                        layout = column(title, a_x, row(p, c))
                        static_layout = column(title, a_x, row(p, c))
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y, c))
                        static_layout = column(title, row(p, a_y, c))
                    else:
                        layout = column(title, row(p, c))
                        static_layout = column(title, row(p, c))

            elif x_track_bed_file:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, j, a_x, row(p, a_y))
                        static_layout = column(title, j, a_x, row(p, a_y))
                    else:
                        layout = column(title, j, a_x, p)
                        static_layout = column(title, j, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, j, row(p, a_y))
                        static_layout = column(title, j, row(p, a_y))
                    else:
                        layout = column(title, j, p)
                        static_layout = column(title, j, p)
            
            else:
                if x_annotation_bed_file:
                    if y_annotation_bed_file:
                        layout = column(title, a_x, row(p, a_y))
                        static_layout = column(title, a_x, row(p, a_y))
                    else:
                        layout = column(title, a_x, p)
                        static_layout = column(title, a_x, p)
                else:
                    if y_annotation_bed_file:
                        layout = column(title, row(p, a_y))
                        static_layout = column(title, row(p, a_y))
                    else:
                        layout = column(title, p)
                        static_layout = column(title, p)

    # Save and show final plot
    html_output_file = output_prefix + ".html"
    output_file(html_output_file)
    
    #show(layout)
    save(layout)
    
    #save layout as png
    png_output_file = output_prefix + ".png"
    export_png(static_layout, filename=png_output_file, scale_factor=3)

    #save layout as svg
    svg_output_file = output_prefix + ".svg"
    export_svg(static_layout, filename=svg_output_file, timeout=1200)


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
    highlighted regions; quantitative bed files wh-
    ch can be loaded as tracks (e.g. read depth); &
    annotation bed files which can be loaded as tr-
    acks (e.g. gene/repeat location bed files).
          
    Please contact ak37@sanger.ac.uk for issues, 
    questions and requests.
    ________________________________________________
                
                  version 2.0.0
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
    parser.add_argument("--paf_file", metavar="FILE", help="minimap2 paf file generated with -c flag or FastGA paf file with the 'de' tag available; x-axis fasta should be the query sequence and y-axis fasta should be the subject sequence (i.e. minimap2/FastGA x.fasta y.fasta)")
    parser.add_argument("--x_index_file", metavar="FILE", help="samtools .fai file of fasta to be on x-axis of dotplot")
    parser.add_argument("--y_index_file", metavar="FILE", help="samtools .fai file of fasta to be on y-axis of dotplot")
    parser.add_argument("--threshold", metavar="FLOAT", type=float, default=90, help="(optional) minimum nucleotide identity of matches to be plotted; default is 90")
    parser.add_argument("--size_threshold", metavar="FLOAT", type=float, default=1000, help="(optional) minimum length of matches to be plotted; default is 1000")
    parser.add_argument("--plot_title", metavar="STR", default="HyraxDotPlot", help="(optional) plot title; default is HyraxDotPlot")
    parser.add_argument("--output_prefix", metavar="STR", default="hyraxdotplot", help="(optional) output file name prefix; default is hyraxdotplot")
    parser.add_argument("--x_feature_bed_file", metavar="FILE", help="(optional) qualitative bed file for x-axis fasta to be loaded as highlighted regions on dot plot")
    parser.add_argument("--y_feature_bed_file", metavar="FILE", help="(optional) qualitative bed file for y-axis fasta to be loaded as highlighted regions on dot plot")
    parser.add_argument("--plot_width", metavar="INT", type=int, default=800, help="(optional) plot width; default is 800")
    parser.add_argument("--plot_height", metavar="INT", type=int, default=600, help="(optional) plot height; default is 600")
    parser.add_argument("--x_track_bed_file", metavar="FILE", help="quantitative bed file to be loaded as a track for x-axis fasta, e.g. windowed coverage bed file")
    parser.add_argument("--x_track_title", metavar="STR", default="Feature Plot", help="(optional) plot title for track loaded on x-axis fasta")
    parser.add_argument("--x_track_feature_name", metavar="STR", default="Feature", help="(optional) feature name for track loaded on x-axis fasta")
    parser.add_argument("--x_track_colour", metavar="STR", default="#56B4E9", help="(optional) colour of track loaded on x-axis fasta; default is #56B4E9")
    parser.add_argument("--y_track_bed_file", metavar="FILE", help="quantitative bed file to be loaded as a track for y-axis fasta, e.g. windowed coverage bed file")
    parser.add_argument("--y_track_title", metavar="STR", default="Feature Plot", help="(optional) plot title for track loaded on y-axis fasta")
    parser.add_argument("--y_track_feature_name", metavar="STR", default="Feature", help="(optional) feature name for track loaded on y-axis fasta")
    parser.add_argument("--y_track_colour", metavar="STR", default="#f49ac2", help="(optional) colour of track loaded on y-axis fasta; default is #f49ac2")
    parser.add_argument("--x_annotation_bed_file", metavar="FILE", help="(optional) gene or repeat annotation bed file for x-axis fasta")
    parser.add_argument("--y_annotation_bed_file", metavar="FILE", help="(optional) gene or repeat annotation bed file for y-axis fasta")
    parser.add_argument("--curation_mode", action='store_true', help="(optional) flag that adds tap tool and allows you to select sequences of interest")

    args = parser.parse_args()

    # Check if positional arguments are missing and print an error message if they are
    
    if not args.coords_file and not args.paf_file:
        print("Error: Please submit a nucmer coords file or a paf file with the de tag available")
    
    missing_args = []
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

    plot_interactive_bokeh(coords_file=args.coords_file, paf_file=args.paf_file, threshold=args.threshold, size_threshold=args.size_threshold,x_cumulative_length_dict=xcumulative_length_dict, y_cumulative_length_dict=ycumulative_length_dict, plot_title=args.plot_title, output_prefix=args.output_prefix, plot_width=args.plot_width, plot_height=args.plot_height, x_feature_bed_file=args.x_feature_bed_file, y_feature_bed_file=args.y_feature_bed_file, x_track_bed_file=args.x_track_bed_file, x_track_title=args.x_track_title, x_track_feature_name=args.x_track_feature_name, x_track_colour=args.x_track_colour, y_track_bed_file=args.y_track_bed_file, y_track_title=args.y_track_title, y_track_feature_name=args.y_track_feature_name, y_track_colour=args.y_track_colour, x_annotation_bed_file=args.x_annotation_bed_file, y_annotation_bed_file=args.y_annotation_bed_file, curation_mode=args.curation_mode)

    
    print("""Plotting finished, enjoy your plot""")

if __name__ == "__main__":
    main()
