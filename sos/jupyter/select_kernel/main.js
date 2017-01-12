

define([
    'jquery',
    'base/js/namespace',
    'base/js/dialog',
    'base/js/utils',
    'services/config'
], function($, Jupyter, dialog, utils, configmod) {


//https://gist.github.com/damianavila/5793132
	(function (IPython) {
	    "use strict";

	    var CellToolbar = IPython.CellToolbar;
	    var slideshow_preset = [];

	    var select_type = CellToolbar.utils.select_ui_generator([
	            ["sos"        ,"sos"        ],
	            ["R"    ,"R"     ],
	            ["python"     ,"python"     ],
	            ],
	            // setter
	            function(cell, value){
	                // we check that the slideshow namespace exist and create it if needed
	                if (cell.metadata.kernel == undefined){cell.metadata.kernel = {}}
	                	cell.metadata.kernel = value
	                	changeStyleOnCellKernel(cell,value)
	                },
	            //geter
	            function(cell){ var ns = cell.metadata.kernel;
	                return (ns == undefined)? undefined: ns.kernel
	                },
	            "");

	    CellToolbar.register_callback('slideshow.select',select_type);

	    slideshow_preset.push('slideshow.select');
	    
	    var reveal_preset = slideshow_preset.slice();
	    

	    CellToolbar.register_preset('Select cell kernel',reveal_preset);
	    console.log('Select cell kernel loaded.');
	    CellToolbar.global_show();
        CellToolbar.activate_preset('Select cell kernel');
         

	}(IPython));


    "use strict";

    // create config object to load parameters
    var base_url = utils.get_body_data("baseUrl");
    var config = new configmod.ConfigSection('notebook', {base_url: base_url});

    config.loaded.then(function() {
        var dropdown = $("<select></select>").attr("id", "kernel_selector")
                                             .css("margin-left", "0.75em")
                                             .attr("class", "form-control select-xs")
                                             // .change(select_kernel);
        Jupyter.toolbar.element.append(dropdown);
        $.each(["sos","R","python"], function(key,value) {   
				     $('#kernel_selector')
				         .append($("<option></option>")
				                    .attr("value",value)
				                    .text(value)); 
				});

        $('#kernel_selector').change(function(){
        	var kernel_type = $("#kernel_selector").val();
	        // var cell = IPython.notebook.get_selected_cell();
	        // var cells = IPython.notebook.get_cells();
	        // var cell = cells[cells.length-1]
	        // cell.metadata.kernel=kernel_type
	        // changeStyleOnCellKernel(cell,kernel_type)
	        // console.log(cell.metadata.kernel)
	        window.default_kernel=kernel_type
    	});


        // var cells = IPython.notebook.get_cells();
        // for (var i = cells.length - 1; i >= 0; --i) {
        //  	 var cell_dropdown = $("<select></select>").attr("id", i+"_cell_kernel_selector")
        //                                      .css("margin-left", "0.75em")
        //                                      .attr("class", "form-control select-xs")
                                           
         	// cells[i].element[0].getElementsByClassName('inner_cell')[0].append(cell_dropdown);
      
         	// Jupyter.notebook._insert_element_at_index(cell_dropdown,i)
         // 	$.each(["sos","R","python"], function(key,value) {   
				     // $('#'+i+"_cell_kernel_selector")
				     //     .append($("<option></option>")
				     //                .attr("value",value)
				     //                .text(value)); 
         // 	});
         // 	$('#'+i+"_cell_kernel_selector").change(function(){
         // 		var kernel_type=$(this).val()
		       //  // var kernel_type = $('#'+i+"_cell_kernel_selector").val();
		       //  console.log(kernel_type)
		       //  // var cell = IPython.notebook.get_selected_cell();
		       //  var cell = Jupyter.notebook.get_cell(i)
		       //  console.log(i)
		       //  cell.metadata.kernel=kernel_type
		       //  changeStyleOnKernel(cell,kernel_type)
		       //  console.log(cell.metadata.kernel)

         // 	})
         // }	

    });

    // will be called when the nbextension is loaded
    function change_kernel() {
        	config.load(); // trigger loading config parameters
        	load_css();
        };
    
    function changeStyleOnCellKernel(cell,type){          
            // this should be loaded from language css file
            cell.element.css('background-color', BC[type]);
            cell.element[0].getElementsByClassName('input_area')[0].style.backgroundColor = BC[type];
            // var original_text=cell.element[0].getElementsByClassName("input_prompt")[0].textContent;
            // original_text=original_text.split(/(\s+)/);
            // console.log(original_text)
            // cell.element[0].getElementsByClassName("input_prompt")[0].textContent=type+" [ ]:";
            // console.log(type+original_text);
            // console.log(cell.element[0].getElementsByClassName("input_prompt")[0].textContent)
        }

     var load_css = function () {
		    var link = document.createElement("link");
		    link.type = "text/css";
		    link.rel = "stylesheet";
		    link.href = require.toUrl("/nbextensions/select_kernel/main.css");
		    document.getElementsByTagName("head")[0].appendChild(link);
	  };

    // return public methods
    return {
        load_ipython_extension : change_kernel
    };
});


// define([
//     'base/js/namespace'
//     ], function(
// 	        Jupyter
// 	) {
// 	    function load_ipython_extension() {

// 	        var handler = function () {
// 				alert('this is an alert from test extension!');
// 			};

// 	        var action = {
// 				    icon: 'fa-comment-o', 
// 				    help    : 'Select kernel',
// 				    help_index : 'zz',
// 				    handler : handler
// 			    };
// 		        var prefix = 'my_extension';
// 			    var action_name = 'select-kernel';
// 				var full_action_name = Jupyter.actions.register(action, name, prefix); 
// 			}

// 		return {
// 				load_ipython_extension: load_ipython_extension
// 		};
//     });
