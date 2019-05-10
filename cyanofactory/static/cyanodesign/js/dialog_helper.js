define(["require", "exports", "jquery", "selectize"], function (require, exports, $) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    class DialogHelper {
        static updateTips(o, txt) {
            $(o).parent(".form-group").addClass("has-error");
            $(o).siblings(".help-block").remove();
            $("<div></div>").addClass("help-block").text(txt).insertAfter(o);
        }
        static checkLength(o, n, min, max) {
            if (o.value.length > max || o.value.length < min) {
                o.classList.add("ui-state-error");
                DialogHelper.updateTips(o, "Length of " + n + " must be between " +
                    min + " and " + max + ".");
                return false;
            }
            else {
                return true;
            }
        }
        static checkBool(o, bool, error) {
            if (bool) {
                return true;
            }
            else {
                DialogHelper.updateTips(o, error);
                return false;
            }
        }
        static checkRegexpPos(o, regexp, n) {
            if (regexp.test(o.value)) {
                DialogHelper.updateTips(o, n);
                return false;
            }
            else {
                return true;
            }
        }
        static checkRegexp(o, regexp, n) {
            if (!(regexp.test(o.value))) {
                DialogHelper.updateTips(o, n);
                return false;
            }
            else {
                return true;
            }
        }
        static checkId(o) {
            let val = o.value;
            if (val.length == 0 || !/[_a-zA-Z]/.test(val[0])) {
                DialogHelper.updateTips(o, "Identifier must start with a letter or '_'");
                return false;
            }
            else if (!/^[_a-zA-Z][_0-9a-zA-Z]*$/.test(val)) {
                DialogHelper.updateTips(o, "Identifier may only contain alphanumeric (a-z, A-Z, 0-9) and '_'");
                return false;
            }
            return true;
        }
        static getFloat(value) {
            if (/^([-+])?([0-9]+(\.[0-9]+)?)$/.test(value))
                return Number(value);
            return NaN;
        }
        ;
    }
    exports.DialogHelper = DialogHelper;
    class ElementWrapper {
        constructor(element) {
            this.html_element = element;
            if (element.tagName == "SELECT") {
                this.type = "text";
            }
            else {
                this.type = element.attributes["type"].value;
            }
            if (this.type != "text" &&
                this.type != "checkbox" &&
                this.type != "radio") {
                throw new Error("Unsupported input type " + this.type);
            }
        }
        static byClass(classname, parent) {
            let elem = parent.getElementsByClassName(classname)[0];
            if (elem == null) {
                throw new Error("ElementWrapper: " + classname + " not found");
            }
            if (elem.tagName != "INPUT" && elem.tagName != "SELECT") {
                throw new Error("ElementWrapper: " + classname + " must be input element but is " + elem.tagName);
            }
            return new ElementWrapper(elem);
        }
        get value() {
            if (this.type == "text") {
                return this.html_element.value;
            }
            else if (this.type == "checkbox" || this.type == "radio") {
                return this.html_element.checked;
            }
            throw new Error("BUG: Unreachable code path taken");
        }
        set value(val) {
            if (this.type == "text") {
                this.html_element.value = val;
            }
            else if (this.type == "checkbox" || this.type == "radio") {
                this.html_element.checked = val;
            }
            else {
                throw new Error("BUG: Unreachable code path taken");
            }
        }
        get element() {
            return this.html_element;
        }
    }
    exports.ElementWrapper = ElementWrapper;
});
