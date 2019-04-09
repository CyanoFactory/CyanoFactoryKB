import * as $ from "jquery";

export class DialogHelper {
    static updateTips(o: JQuery<HTMLElement>, txt: string) {
            o.parent(".form-group").addClass("has-error");
            o.siblings(".help-block").remove();
            $("<div></div>").addClass("help-block").text(txt).insertAfter(o);
    }

    static checkLength(o: JQuery<HTMLElement>, n: string, min: number, max: number) {
        if (o.val().length > max || o.val().length < min) {
            o.addClass("ui-state-error");
            DialogHelper.updateTips(o, "Length of " + n + " must be between " +
                min + " and " + max + ".");
            return false;
        } else {
            return true;
        }
    }

    static checkBool(o: JQuery<HTMLElement>, bool: boolean, error: string) {
        if (bool) {
            return true;
        } else {
            DialogHelper.updateTips(o, error);
            return false;
        }
    }

    static checkRegexpPos(o: JQuery<HTMLElement>, regexp: RegExp, n) {
        if ( regexp.test(<string>o.val()) ) {
            DialogHelper.updateTips(o, n);
            return false;
        } else {
            return true;
        }
    }

    static checkRegexp(o: JQuery<HTMLElement>, regexp: RegExp, n) {
        if (!( regexp.test(<string>o.val()) )) {
            DialogHelper.updateTips(o, n);
            return false;
        } else {
            return true;
        }
    }
}