waitIndicator = (function () {
    var pleaseWaitDiv = $('<div class="please-wait" style="position: absolute;top:0;left:0;width: 100%;height:100%;z-index:2;opacity:0.4;background-color:grey"></div>');
    return {
        show: function(target) {
            var content = $(target);
            var spinner = new Spinner().spin();
            content.append(pleaseWaitDiv);
            content.append(spinner.el);
        },
        hide: function (target) {
            var content = $(target);
            content.find(".spinner").remove();
            content.find(".please-wait").remove();
        }
    };
})();
