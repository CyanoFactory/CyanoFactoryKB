$(document).ready(function() {
    $("form#search_form").submit(function(e) {
        var search_box = $("form textarea#search_box")
        var text = search_box.val()
        text = text.replaceAll("#", "%23")
        search_box.val(text)
    });

    $("form#load_form").submit(function(e) {
        e.preventDefault();
    });

    $("form#save_form").submit(function(e) {
        e.preventDefault();
        $.ajax({
            url: "{% url "kegg.views.index_ajax" %}",
            context: document.body,
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

    $("form#load_form input#load_button").click(function(e) {
        $.ajax({
            url: "{% url "kegg.views.index_ajax" %}",
            context: document.body,
            data: {"op": "load", "pk": $("form#load_form select").children(":selected").attr("name")}
        }).done(function(result) {
            var json = jQuery.parseJSON(result);
            var search_box = $("form#search_form textarea#search_box");
            search_box.val(json.query);
            $("form#save_form input#save_input").val(json.name);
        });
    });

    $("form#load_form input#delete_button").click(function(e) {
        $.ajax({
            url: "{% url "kegg.views.index_ajax" %}",
            context: document.body,
            data: {"op": "delete", "pk": $("form#load_form select").children(":selected").attr("name")}
        }).done(function(result) {
            $("form#load_form select").children(":selected").remove();
        });
    });
});
