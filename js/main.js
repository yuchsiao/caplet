$(document).ready(function(){
    $('.bxslider').bxSlider({
        mode: 'fade',
        captions: true
    });

    // $('footer').append($(document).height);
//    $(document).keydown(function(e){
//        if (e.keyCode == 37) {
//            alert( "left pressed" );
//            slider.goToPrevSlide();
//            return false;
//        }
//    });

    var set_height = function(){
        var content_height = $('#banner').height()
                            +$('#main').height()
                            +$('#footer').height() + 100;
        // alert('content = ' + content_height + '\nwindow = ' + $(window).height());
        if ( content_height < $(window).height() ){
            $('#center-column').height($(window).height());
        }
        else{
            $('#center-column').height(content_height); 
        }
        $('#main').show();
    };

    var set_height_given = function(content_height){
        if ( content_height < $(window).height() ){
            $('#center-column').height($(window).height());
        }
        else{
            $('#center-column').height(content_height); 
        }
    };

    var load_content = function(id, filename){
        return $.get(filename, function(data){
            //* Publish md file into html
            $(id).html(markdown.toHTML(data));
            //* Add attr to link to open in new windows
            $(id+' a').attr('target', '_blank');
        }, 'text');
    };


    var tab_color_current = "#E65CB8";
    var tab_color_other   = "#00b7d1";

    var tabs = new Array();
    tabs[0] = "index";
    tabs[1] = "news";
    tabs[2] = "feature";
    tabs[3] = "screenshot";
    tabs[4] = "download";
    tabs[5] = "installation";
    tabs[6] = "tutorial";
    tabs[7] = "publication";

    var load_main = function(url){
        //* Color tab buttons
        for (var i=0; i<tabs.length; ++i){
            if (tabs[i]+'.html'==url){
                $('#'+tabs[i]+'-tab').css('background-color', tab_color_current);
                // alert('')
            }
            else{
                $('#'+tabs[i]+'-tab').css('background-color', tab_color_other);
            }
        }

        //* Load pages
        $('#main').load(url + ' #main', function(){
            if (url != 'screenshot.html'){
                $('#main').hide();
            }

            //* ScreenShot content
            if ( url == 'screenshot.html' ){
                var slider = $('.bxslider').bxSlider({
                    mode: 'fade',
                    captions: true,
                    speed: 100
                });
                $(document).keydown(function(e){
                    if (e.keyCode == 37) {
                        slider.goToPrevSlide();
                        return false;
                    }
                    else if (e.keyCode == 39) {
                        slider.goToNextSlide();
                        return false;
                    }
                });
                set_height_given(1000);
            }

            //* News content
            else if ( url == 'news.html'){
                $.when( load_content('#news-markdown'   , 'markdown/NEWS.md') , 
                        load_content('#release-markdown', 'markdown/RELEASE.md') )
                            .done(function(){
                    set_height();
                });
            }

            //* Installation content
            else if ( url == 'installation.html'){
                $.when( load_content('#requirement-markdown',  'markdown/REQUIREMENT.md'),
                        load_content('#installation-markdown', 'markdown/INSTALLATION.md') , 
                        load_content('#testenv-markdown', 'markdown/TESTENV.md') )
                            .done(function(){
                    set_height();
                });
            }

            //* Tutorial content
            else if ( url == 'tutorial.html'){
                $.when( load_content('#quickstart-markdown', 'markdown/QUICKSTART.md'),
                        load_content('#limitation-markdown', 'markdown/LIMITATION.md'),
                        load_content('#format-markdown', 'markdown/FORMAT.md') )
                            .done(function(){
                    set_height();
                });
            }

            //* Feature
            else if ( url == 'feature.html' ){
                set_height();
            }

            //* Publication
            else if ( url == 'publication.html' ){
                set_height();
            }

            //* Download
            else if ( url == 'download.html' ){
                set_height();
            }

            //* The rest
            else{
                set_height();
            }

        });
    };


    set_height();

    $(window).resize(set_height);

    $('#main, nav').on('click', 'a.tab-link', function(){
        var url=$(this).attr('href');
        load_main(url);
        return false;
    });

    for (var i=0; i<tabs.length; ++i ){
        if($('#'+tabs[i]).length){
            load_main(tabs[i]+'.html');
            break;
        }        
    }
    // if($('#news').length){
    //     load_main('news.html');        
    // }
    // else if($('#feature').length){
    //     load_main('feature.html');        
    //     // alert('here');
    // }
    // else if($('#screenshot').length){
    //     load_main('screenshot.html');        
    // }
    // else if($('#tutorial').length){
    //     load_main('tutorial.html');        
    // }

}); // end ready


//$('.bxslider').bxSlider({
//    mode: 'fade',
//    captions: true,
//    infiniteLoop: false,
//    hideControlOnEnd: true
//});
