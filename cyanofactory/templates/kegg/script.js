$(document).ready(function() {
    $("form#load_form").submit(function(e) {
        e.preventDefault();
    });

    $("form#save_form").submit(function(e) {
        e.preventDefault();
        $.ajax({
            url: '{% url "kegg.views.index_ajax" %}',
            context: document.body,
            method: "POST",
            data: {
                "op": "save",
                "name": $("form#save_form input#save_input").val(),
                "query": $("form#search_form textarea#search_box").val()
            }
        }).done(function(result) {
            var json = jQuery.parseJSON(result);
            var load_item = $("form#load_form select").children("select");
            if (json.created) {
                load_item.end().append('<option selected="selected" name="' + json.pk + '">' + json.name + '</option>');
            }
        });
    });

    $("form#load_form #load_button").click(function(e) {
        $.ajax({
            url: '{% url "kegg.views.index_ajax" %}',
            context: document.body,
            method: "POST",
            data: {"op": "load", "pk": $("form#load_form select").children(":selected").attr("name")}
        }).done(function(result) {
            var json = jQuery.parseJSON(result);
            var search_box = $("form#search_form textarea#search_box");
            search_box.val(json.query);
            $("form#save_form input#save_input").val(json.name);
        });
    });

    $("form#load_form #delete_button").click(function(e) {
        $.ajax({
            url: '{% url "kegg.views.index_ajax" %}',
            context: document.body,
            data: {"op": "delete", "pk": $("form#load_form select").children(":selected").attr("name")}
        }).done(function(result) {
            $("form#load_form select").children(":selected").remove();
        });
    });

    var isDragging = false;

    $("#content").on("mousedown", "a", function(event) {
        isDragging = false;
        $(this).data('page', {x: event.pageX, y: event.pageY})
    });
    $("#content").on("mousemove", "a", function(event) {
        var p = $(this).data('page');
        if (p !== undefined) {
            if (Math.abs(p.x - event.pageX) > 4 ||
                Math.abs(p.y - event.pageY) > 4) {
                isDragging = true;
            }
        }
    });
    $("#content").on("click", "a", function(event) {
        var href = undefined;

        if (typeof $(this).attr("href") !== "undefined") {
            href = $(this).attr("href");
        }
        else if (typeof($(this).attr("xlink:href") !== "undefined")) {
            href = $(this).attr("xlink:href");
        }

        if (href === undefined) {
            return;
        }

        if (href.indexOf('{% url "kegg.views.index" %}') == 0) {
            event.preventDefault();

            if (!isDragging) {
                // redirect search form post target
                var search_form = $("form#search_form");
                search_form.attr("action", href);
                search_form.submit();
            }
        }
    });
});
