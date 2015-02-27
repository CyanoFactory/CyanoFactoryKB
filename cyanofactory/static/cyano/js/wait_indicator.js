waitIndicator = (function () {
    var pleaseWaitDiv = $('<div class="please-wait" style="position: absolute;top:0;left:0;width: 100%;height:100%;z-index:2;opacity:0.4;background-color:grey"></div>');
    return {
        show: function(target) {
            var content = target.find(".modal-content");
            var spinner = new Spinner().spin();
            content.append(pleaseWaitDiv);
            content.append(spinner.el);
            target.on('hide.bs.modal', function(e) {
                e.preventDefault();
            });
        },
        hide: function (target) {
            var content = $(target).find(".modal-content");
            content.find(".spinner").remove();
            content.find(".please-wait").remove();
            target.off('hide.bs.modal');
        }
    };
})();
