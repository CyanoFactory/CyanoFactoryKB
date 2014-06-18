'use strict';

(function($) {
    $(document).ready(function(){
        var page_loading = false;

        $(document).on("click", "a.endless_more", function() {
            if (page_loading) {
                return false;
            }
            page_loading = true;
            var container = $(this).closest(".endless_container");
            var loading = container.find(".endless_loading");
            $(this).hide();
            loading.show();
            var data = "querystring_key=" + $(this).attr("rel").split(" ")[0];
            $.get($(this).attr("href"), data, function(data) {
                container.before(data);
                container.remove();
                page_loading = false;
            });
            return false;
        });
        $(document).on("click", "a.endless_page_link", function() {
            var page_template = $(this).closest(".endless_page_template");
            if (!page_template.hasClass("endless_page_skip")) {
                var data = "querystring_key=" + $(this).attr("rel").split(" ")[0];
                page_template.load($(this).attr("href"), data);
                return false;
            }
            return true;
        });
    });
})(jQuery);
