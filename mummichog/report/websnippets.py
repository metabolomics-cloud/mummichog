# specific for mummichog;
# js libraries linked from web
# last update: 2017-07-09, Shuzhao Li


HTML_HEAD = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"> 
        <html xmlns="http://www.w3.org/1999/xhtml"> 
        <head> 
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" /> 
        <title>Mummichog Report</title> 

        <style type="text/css"> 
        body {width: 960px; padding: 10px; }
        
        div.colorbar {margin-bottom: 7px}
        div.networkvisual { margin: auto; height: 600px; width: 960px; overflow: auto !important;}
        div.svg2 { margin: auto; height: 600px; width: 960px; overflow: auto !important;}

        div.stats {margin-top: 10; font-size:0.7em; color:#888;}

        div.moduleline{font-size: 0.9em; font-weight: bold;
        }
        div.metabolites {
            font-size: 0.7em;
            padding-left: 15px; padding-bottom: 7px;
        }
        
        footer { margin-top: 30px;    padding: 10px 0;    border-top: 1px solid #582E2E;    font-size:0.7em;    color:#888;}
        
        h1, h2, h3, h4, h5, h6 {
            font-family: 'Trebuchet MS', 'Lucida Grande', Arial, Sans-Serif;
            font-weight: bold;
        }
        
        h1 { font-size: 1.4em; }
        h2 { font-size: 1.2em; }
        h3 { font-size: 1em; }
        h4 { font-size: 0.9em;     padding-left: 15px; padding-bottom: 2px; margin-bottom:0px;}
        th {
            background: #DAFDCF;
            font-size:0.9em;
            text-align: left;
            padding:5px;
        }
        
        tr:nth-child(even) {background-color: #f2f2f2}
        tr:hover {background-color: #f5f5f5}
        
        td {
            font-family: Verdana, Arial, sans-serif;
            color: black;
            font-size:0.7em;
            margin:10px;
            margin-top:0px;
            padding: 5px;
        }
        
        .node {
          stroke: #fff;
          stroke-width: 1.5px;
        }
        
        .link {
          stroke: #999;
          stroke-opacity: .6;
        }
        div.network_selection {float: left; margin-top: 6; margin-bottom: 10; font-size:0.8em;}
        </style> 
        <script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>
        <script src="http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js" charset="utf-8"></script>
        <script src="http://mummichog.org/download/cytoscape.min.js" charset="utf-8"></script>
        <script src="https://cdn.rawgit.com/cytoscape/cytoscape.js-spread/1.0.9/cytoscape-spread.js"></script>        

        </head> 
        <body> 
        """
        
HTML_END = """</body> </html> """
        
javascript_HEAD = """
        <script type="text/javascript" charset="utf-8">
        """
        
javascript_END = """
        var w = 960, h = 600;
        var color = d3.scale.category20b();
        var node_sizes = [16, 8, 32]
        
        var scalebar = d3.select("#colorbar").append("svg").attr("id", "svg_colorbar").attr("width", w).attr("height", 10);
        var coloridx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19];
        scalebar.selectAll("rect")
           .data(coloridx)
           .enter()
           .append("rect")
           .attr("x", function(d, i) {return i * 48;})
           .attr("y", 0)
           .attr("width", 48)
           .attr("height", 10)
           .attr("fill", function(d, i) { return color(i); });
           
        
        // initial draw
        cyto_draw_figure(cytonodes[0], cytoedges[0], node_sizes[0]);
    	var cy = $('#networkvisual').cytoscape('get');
    	cy.panBy({ x: 0, y: -200 });
    	
        // support cytoscape.js  
		function cyto_draw_figure(nodes, links, node_size) {  

                    var svg = document.createElement('div');
                    svg.setAttribute("id","svg2");
                    svg.setAttribute("class","svg2");
                    var nw = document.getElementById("networkvisual");
                    nw.appendChild(svg);

				$('#svg2').cytoscape({
						style: cytoscape.stylesheet()
						.selector('node')
						  .css({
							'content': 'data(id)',
							'font-size': node_size,
							'min-zoomed-font-size': 8,
							'background-color': 'mapData(group, 0, 19, blue, red)'
						  })
						.selector('edge')
						  .css({
							'target-arrow-shape': 'triangle',
							'width': 4,
							'line-color': '#ddd',
							'target-arrow-color': '#ddd'
						  })
						.selector('.highlighted')
						  .css({
							'background-color': '#61bffc',
							'line-color': '#61bffc',
							'target-arrow-color': '#61bffc',
							'transition-property': 'background-color, line-color, target-arrow-color',
							'transition-duration': '0.5s'
						  }),
  
					  layout: {
						name: 'spread',
						minDist: 40,
						directed: true,
						padding: 10,
						//fit: false,
					  },
					  
					  elements: {
						  nodes: nodes, 
						  edges: links
						},
  
					  zoom: 1,
					  minZoom: 0.1,
 					  maxZoom: 5,

				});
	
		}
	

        // support d3.js  
        function d3_draw_figure(nodes, links, node_size) {

            var svg = document.createElement('div');
            svg.setAttribute("id","svg2");
            svg.setAttribute("class","svg2");
            var nw = document.getElementById("networkvisual");
            nw.appendChild(svg);

            var svg = d3.select("#svg2").append("svg").attr("id", "svg_networkvisual").attr("width", w).attr("height", h);
            var force = d3.layout.force()
                      .charge(-480)
                      .linkDistance(80)
                      .gravity(0.3)
                      .size([w, h]);
    
            force.nodes(nodes)
              .links(links)
              .start();
            
            var link = svg.selectAll(".link")
              .data(links)
              .enter().append("line")
              .attr("class", "link")
              .style("stroke-width", 2);
            
            var node = svg.selectAll(".node")
              .data(nodes)
              .enter().append("circle")
              .attr("class", "node")
              .attr("r", node_size)
              .style("fill", function(d) { return color(d.group); })
              .style("fill-opacity", .6)
              .call(force.drag);
            
            var text = svg.selectAll(".text")
              .data(nodes)
              .enter().append("text")
              .attr("font-family", "sans-serif")
              .attr("font-size", "9px")
              .attr("text-anchor", "middle")
              .text(function(d) { return d.name; });
            
            force.on("tick", function() {
                link.attr("x1", function(d) { return d.source.x; })
                    .attr("y1", function(d) { return d.source.y; })
                    .attr("x2", function(d) { return d.target.x; })
                    .attr("y2", function(d) { return d.target.y; });
            
                node.attr("cx", function(d) { return d.x; })
                    .attr("cy", function(d) { return d.y; });
            
                text.attr("x", function(d) { return d.x; })
                    .attr("y", function(d) { return d.y; });
              });
        }
    
    
        function inv(i) {
           if (i==0)
              return 0
           else if (i==1)
              return 2
           else
              return 1
}

    	// added support of two styles
        function updateDrawing() {
            var idx = document.getElementById("nmodule_select").selectedIndex;
            var ndx = document.getElementById("nsize_select").selectedIndex;
            var vdx = document.getElementById("vstyle_select").selectedIndex;
            
            //d3.select("#svg_networkvisual").remove();
            var svg = document.getElementById("svg2");
	    svg.parentNode.removeChild(svg);

            //var cy = $('#networkvisual').cytoscape('get');
            //cy.remove( cy.elements() );
            
            // 0 for cytoscape style, drag; 1 for d3 style, force
            if (vdx == 0) {
            cyto_draw_figure(cytonodes[idx], cytoedges[idx], node_sizes[inv(ndx)]);
            
    		var cy = $('#networkvisual').cytoscape('get');
    		cy.panBy({ x: 0, y: -200 });
            
            } else {
            d3_draw_figure(nodes[idx], links[idx], node_sizes[ndx]);
            }
        }
    
    </script>
    """



