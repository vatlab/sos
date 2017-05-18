{%- extends 'full.tpl' -%}


{%- block header -%}
{{ super() }}

 <link rel="stylesheet" href="http://code.jquery.com/ui/1.11.4/themes/smoothness/jquery-ui.css">

<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js"></script>
<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jqueryui/1.9.1/jquery-ui.min.js"></script>

<style>  /* defined here in case the main.css below cannot be loaded */
.lev1 {margin-left: 80px}
.lev2 {margin-left: 100px}
.lev3 {margin-left: 120px}
.lev4 {margin-left: 140px}
.lev5 {margin-left: 160px}
.lev6 {margin-left: 180px}
</style>

<link rel="stylesheet" type="text/css" href="../../css/jt.css">
<link rel="stylesheet" type="text/css" href="../../css/toc2.css">

<script src="../../js/doc_toc.js"></script>

 <script src="../../js/docs.js"></script>

<script>
    MathJax.Hub.Config({
        "HTML-CSS": {
            preferredFont: "TeX",
            availableFonts: ["STIX","TeX"],
            styles: {
                scale: 110,
                ".MathJax_Display": {
                    "font-size": "110%",
                }
            }
        }
    });
</script>

<script>
$( document ).ready(function(){

            var cfg={'threshold':{{ nb.get('metadata', {}).get('toc', {}).get('threshold', '3') }},     // depth of toc (number of levels)
             // 'number_sections': {{ 'true' if nb.get('metadata', {}).get('toc', {}).get('number_sections', False) else 'false' }},  // sections numbering
             'number_sections': false, 
             'toc_cell': false,          // useless here
             'toc_window_display': true, // display the toc window
             "toc_section_display": "block", // display toc contents in the window
             // 'sideBar':{{ 'true' if nb.get('metadata', {}).get('toc', {}).get('sideBar', False) else 'false' }},      
             'sideBar':true,       // sidebar or floating window
             'navigate_menu':false       // navigation menu (only in liveNotebook -- do not change)
            }

            var st={};                  // some variables used in the script
            st.rendering_toc_cell = false;
            st.config_loaded = false;
            st.extension_initialized=false;
            st.nbcontainer_marginleft = $('#notebook-container').css('margin-left')
            st.nbcontainer_marginright = $('#notebook-container').css('margin-right')
            st.nbcontainer_width = $('#notebook-container').css('width')
            st.oldTocHeight = undefined
            st.cell_toc = undefined;
            st.toc_index=0;

            // fire the main function with these parameters



            table_of_contents(cfg,st);

            var file=documentationDict[$("h1:first").attr("id")];
            var path="http://vatlab.github.io/SOS"
            // var path="file:///Users/jma7/Development/SOS/docs"
            $("#toc-level0 a").css("color","#126dce");
            $('a[href="#'+$("h1:first").attr("id")+'"]').hide()
            var docs=documentation;
            var pos=documentation.indexOf(file);
        
            for (var a=pos;a>=0;a--){
                  var name=docs[a]
                  $('<li><a href="'+path+'/doc/documentation/'+name+'.html">'+name.replace(/_/g," ")+'</a></li>').insertBefore("#toc-level0 li:eq(0)");
            }
            $('a[href="'+path+'/doc/documentation/'+file+'.html'+'"]').css("color","#126dce");


            // $('<li id="indexHome"><a href="/Users/jma7/Development/SOS/docs/index.html#documentation"><b>Home<b></a></li>').insertBefore("#toc-level0 li:eq(0)");
            for (var a=pos+1;a<docs.length;a++){
                  var name=docs[a]
                  $(".toc #toc-level0").append('<li><a href="'+path+'/doc/documentation/'+name+'.html">'+name.replace(/_/g," ")+'</a></li>');
            }

            // var path="file:///Users/jma7/Development/SOS/website"
            // $(".toc #toc-level0").append('<li id="indexHome"><a href="'+path+'/index.html" ><b>Home<b></a></li>');

            // var docs=documentation
            // for (var a =0;a<docs.length;a++){
            //       var name =docs[a];
            //       $(".toc #toc-level0").append('<li><a href="'+path+'/doc/documentation/'+name+'.html">'+name.split("_").join(" ")+'</a></li>');
            // }
            // var home=$("#toc-level0 #indexHome");
          
            // home.insertBefore("#toc-level0 li:eq(0)");

            // $("#toc-level0 li").filter(".home").insertBefore($("#toc-level0 li").filter(':nth-child(1)'));
            // $("#toc").attr("style","max-height:938px")


    });
</script>


{%- if nb['metadata'].get('sos',{}).get('kernels',none) is not none -%}

<style>  /* defined here in case the main.css below cannot be loaded */

 {% for item in nb['metadata'].get('sos',{}).get('kernels',{}) %}

{%- if item[2] -%}
.lan_{{item[0]}} {background-color: {{item[2]}} !important }

{%- else -%}
.lan_{{item[0]}} {}

{%- endif -%}

{% endfor %}

</style>


{%- endif -%}

{%- endblock header -%}


{% block codecell %}

{%- if cell['metadata'].get('kernel',none) is not none -%}
	<div class="cell border-box-sizing code_cell rendered lan_{{cell['metadata'].get('kernel', none)}}">
	{{ super() }}
	</div>
{%- else -%}
	{{ super() }}
{%- endif -%}

{%- endblock codecell %}
