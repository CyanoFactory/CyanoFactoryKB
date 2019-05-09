import * as $ from "jquery";

export class DialogHelper {
    static updateTips(o: HTMLElement, txt: string) {
            $(o).parent(".form-group").addClass("has-error");
            $(o).siblings(".help-block").remove();
            $("<div></div>").addClass("help-block").text(txt).insertAfter(o);
    }

    static checkLength(o: HTMLElement, n: string, min: number, max: number) {
        if ((<any>o).value.length > max || (<any>o).value.length < min) {
            o.classList.add("ui-state-error");
            DialogHelper.updateTips(o, "Length of " + n + " must be between " +
                min + " and " + max + ".");
            return false;
        } else {
            return true;
        }
    }

    static checkBool(o: HTMLElement, bool: boolean, error: string) {
        if (bool) {
            return true;
        } else {
            DialogHelper.updateTips(o, error);
            return false;
        }
    }

    static checkRegexpPos(o: HTMLElement, regexp: RegExp, n) {
        if ( regexp.test(<string>(<any>o).value)) {
            DialogHelper.updateTips(o, n);
            return false;
        } else {
            return true;
        }
    }

    static checkRegexp(o: HTMLElement, regexp: RegExp, n) {
        if (!( regexp.test(<string>(<any>o).value))) {
            DialogHelper.updateTips(o, n);
            return false;
        } else {
            return true;
        }
    }

    static checkId(o: HTMLElement) {
        let val: string = (<any>o).value;
        if (val.length == 0 || !/[_a-zA-Z]/.test(val[0])) {
            DialogHelper.updateTips(o, "Identifier must start with a letter or '_'");
            return false;
        } else if (!/^[_a-zA-Z][_0-9a-zA-Z]*$/.test(val)) {
            DialogHelper.updateTips(o, "Identifier may only contain alphanumeric (a-z, A-Z, 0-9) and '_'");
            return false;
        }
        return true;
    }

    static getFloat(value: string): float {
        if (/^([-+])?([0-9]+(\.[0-9]+)?)$/.test(value))
            return Number(value);
        return NaN;
    };
}

export class ElementWrapper<T> {
    private readonly html_element: HTMLInputElement;
    private readonly type: string;

    static byClass<T>(classname: string, parent: HTMLElement): ElementWrapper<T> {
        let elem: HTMLElement | null = <HTMLElement>parent.getElementsByClassName(classname)[0];
        if (elem == null) {
            throw new Error("ElementWrapper: " + classname + " not found");
        }

        if (elem.tagName != "INPUT") {
            throw new Error("ElementWrapper: " + classname + " must be input element")
        }

        return new ElementWrapper<T>(<HTMLInputElement>elem);
    }

    constructor(element: HTMLInputElement) {
        this.html_element = element;
        this.type = element.attributes["type"].value;

        if (this.type != "text" &&
            this.type != "checkbox") {
            throw new Error("Unsupported input type " + this.type);
        }
    }

    get value(): T {
        if (this.type == "text") {
            return <T>(<any>this.html_element).value;
        } else if (this.type == "checkbox") {
            return <T>(<any>this.html_element).checked;
        }

        throw new Error("BUG: Unreachable code path taken");
    }

    set value(val: T) {
        if (this.type == "text") {
            (<any>this.html_element).value = val;
        } else if (this.type == "checkbox") {
            (<any>this.html_element).checked = val;
        } else {
            throw new Error("BUG: Unreachable code path taken");
        }
    }

    get element(): HTMLElement {
        return this.html_element;
    }
}
