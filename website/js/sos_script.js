$(document).ready(function(){
	$("#intro_content").load("https://raw.githubusercontent.com/BoPeng/SOS/onepage/doc/README.html"); 
	
	// var tutorials=["Quick_Start_Guide","SoS-QuickStart",
	// 				"SoS_Docker_Guide","Using_SoS_with_IPython","Using_SoS_with_Spyder"]
	// $("#tutorial > .container").append('<div class="row">')
	// for (var a =0;a<tutorials.length;a++){
	// 	var name =tutorials[a];
	// 	var oneString='<div class="col-md-4 col-sm-6 portfolio-item">'
 //        				+'<div class="portfolio-caption">';
	// 	oneString+='<a href="./doc/tutorial_html/'+name+'.html" class="portfolio-link"><h4>'+name+'</h4></a>';
	// 	oneString+='</div></div>';       	
	// 	 $("#tutorial > .container").append(
	//        	oneString
	//      );
	// }	
	// $("#tutorial > .container").append('</div>');	



	// var documentations=["Auxiliary_Steps","Command_Line_Options","Configuration_Files","Extending_SoS","File_Structure",
	// 				"Overview_of_SoS","SoS_Actions","SoS_Syntax","Step_Dependencies","Step_Input","Step_Options",
	// 				"Step_Output","Step_Process_and_External_Tasks","StringInterpolation","UnderstandingTargets",
	// 				"User_Interface","WorkflowSpecification"]
	// $("#documentation > .container").append('<div class="row">')
	// for (var a =0;a<documentations.length;a++){
	// 	var name =documentations[a];
	// 	var oneString='<div class="col-md-4 col-sm-6 portfolio-item">'
 //        				+'<div class="portfolio-caption">';
	// 	oneString+='<a href="./doc/documentation_html/'+name+'.html" class="portfolio-link"><h4>'+name+'</h4></a>';
	// 	oneString+='</div></div>';       	
	// 	 $("#documentation > .container").append(
	//        	oneString
	//      );
	// }	
	// $("#documentation > .container").append('</div>');		








	var dir = "../doc/documentation_html/";
	var fileextension = ".html";
	$("#documentation > .container").append('<div class="row">')
	$.ajax({
	    //This will retrieve the contents of the folder if the folder is configured as 'browsable'
	    url: dir,
	    success: function (data) {
	        //List all .png file names in the page
	        $(data).find("a:contains(" + fileextension + ")").each(function () {
	
	        	var cons=this.href.split("/");
	        	var name=cons.slice(1).slice(-2);
	        	var path="./"+name[0]+"/"+name[1];
	            var oneString='<div class="col-md-4 col-sm-6 portfolio-item">'
        				+'<div class="portfolio-caption">';
        		oneString+="<a href="+path+">"+name[1].replace("html","")+"</a>";
        		oneString+='</div></div>';       	
        		 $("#documentation > .container").append(
			       	oneString
			     );
	  
	        });
	    }
	});
	$("#documentation > .container").append('</div>');

	var dir = "../doc/tutorial_html/";
	var fileextension = ".html";
	$("#tutorial > .container").append('<div class="row">')
	$.ajax({
	    //This will retrieve the contents of the folder if the folder is configured as 'browsable'
	    url: dir,
	    success: function (data) {
	        //List all .png file names in the page
	        $(data).find("a:contains(" + fileextension + ")").each(function () {
	
	        	var cons=this.href.split("/");
	        	var name=cons.slice(1).slice(-2);
	        	var path="./"+name[0]+"/"+name[1];
	            var oneString='<div class="col-md-4 col-sm-6 portfolio-item">'
        				+'<div class="portfolio-caption">';
        		oneString+="<a href="+path+">"+name[1].replace("html","")+"</a>";
        		oneString+='</div></div>';       	
        		 $("#tutorial > .container").append(
			       	oneString
			     );
	  
	        });
	    }
	});
	$("#tutorial > .container").append('</div>');
	
});
