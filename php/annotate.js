//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

var http_req;

function toggle_annotation(id) {
    var div_obj = document.getElementById(id + "_div"); 
    var a_obj   = document.getElementById(id + "_ann"); 
    var frm_obj = document.getElementById(id + "_frm"); 

    // Turn ON the anchor, turn OFF the form
    if (a_obj.style.display == "none") {
        a_obj.style.display = "";
        div_obj.style.display = "none";
        if (a_obj.innerHTML.length > 0) 
            a_obj.innerHTML = "annotate";
    }
    // Turn OFF the anchor, turn ON the form
    else {
        a_obj.style.display   = "none";
        div_obj.style.display = "";
        if (a_obj.innerHTML != "annotate") 
            frm_obj.ext_id.value  = a_obj.innerHTML;
    }
}

function annotate_marker(id) {
    //
    // Fetch the marker annotation
    //
    var form_obj = document.getElementById(id + "_frm"); 
    var url = 
        form_obj.url.value + "&" + 
        "ext_id=" + form_obj.ext_id.value;

    // 
    // Prepare and send XMLHttpRequest Object.
    //
    http_req = false;

    try {
	http_req = new XMLHttpRequest();
    } catch(e) {
	http_req = false;
    }

    if (http_req) {
	http_req.onreadystatechange = process_annotation;
	http_req.open("GET", url, true);
	http_req.send("");
    }

    toggle_annotation(id);
}

function process_annotation() {
    //
    // Possible readyState values:
    // 0 = uninitialized
    // 1 = loading
    // 2 = loaded
    // 3 = interactive
    // 4 = complete
    //
    if (http_req.readyState == 4) {

        // Check that the status is "OK"
        if (http_req.status == 200) {

            var xml_doc = http_req.responseXML;
            var tag_obj = xml_doc.getElementsByTagName("marker_id");
            var txt_obj = xml_doc.getElementsByTagName("text");
            var tag_id  = tag_obj[0].childNodes[0].nodeValue;
            var txt;
            if (txt_obj[0].childNodes.length > 0)
                txt = txt_obj[0].childNodes[0].nodeValue;
            else
                txt = "annotate";

            var a_obj = document.getElementById(tag_id + "_ann");
            a_obj.innerHTML = txt;

        } else {
            alert("There was a problem retrieving the XML data:\n" + http_req.statusText);
        }
    }
}

function toggle_correction(id) {
    var div_obj = document.getElementById(id + "_div");
    var sel_obj = document.getElementById(id + "_sel");
    var s_obj   = document.getElementById(id);

    if (div_obj.style.display == "none") {
        div_obj.style.display = "";
	sel_obj.style.display = "none";
    }
    else {
        div_obj.style.display = "none";
	sel_obj.style.display = "";
	s_obj.focus();
    }
}

function cancel_correction(id) {
    var div_obj = document.getElementById(id + "_div");
    var sel_obj = document.getElementById(id + "_sel");

    div_obj.style.display = "";
    sel_obj.style.display = "none";
}

function correct_genotype(id, url) {
    //
    // Fetch the marker annotation
    //
    var sel_obj = document.getElementById(id); 
    url = url + "&" + "gtype=" + sel_obj.options[sel_obj.selectedIndex].text;

    // 
    // Prepare and send XMLHttpRequest Object.
    //
    http_req = false;

    try {
	http_req = new XMLHttpRequest();
    } catch(e) {
	http_req = false;
    }

    if (http_req) {
	http_req.onreadystatechange = process_correction;
	http_req.open("GET", url, true);
	http_req.send("");
    }

    toggle_correction(id);
}

function process_correction() {
    //
    // Possible readyState values:
    // 0 = uninitialized
    // 1 = loading
    // 2 = loaded
    // 3 = interactive
    // 4 = complete
    //
    if (http_req.readyState == 4) {

        // Check that the status is "OK"
        if (http_req.status == 200) {

            var xml_doc = http_req.responseXML;
            var tag_obj = xml_doc.getElementsByTagName("div_id");
            var div_id  = tag_obj[0].childNodes[0].nodeValue;
            var cor_obj = xml_doc.getElementsByTagName("corrected");
            var cor     = cor_obj[0].childNodes[0].nodeValue;
            var txt_obj = xml_doc.getElementsByTagName("gtype");
            var gtype   = txt_obj[0].childNodes[0].nodeValue;

            var div_obj = document.getElementById(div_id + "_div");
	    var txt;
	    if (cor == "true")
		txt = "<span class=\"corrected\">" + gtype + "</span>";
	    else
		txt = gtype;

            div_obj.innerHTML = "<a onclick=\"toggle_correction('" + div_id + "')\">" + txt + "</a>";

        } else {
            alert("There was a problem retrieving the XML data:\n" + http_req.statusText);
        }
    }
}

function toggle_population(id) {
    var div_obj = document.getElementById(id + "_div"); 
    var a_obj   = document.getElementById(id + "_pop"); 
    var frm_obj = document.getElementById(id + "_frm"); 

    // Turn ON the anchor, turn OFF the form
    if (a_obj.style.display == "none") {
        a_obj.style.display = "";
        div_obj.style.display = "none";
    }
    // Turn OFF the anchor, turn ON the form
    else {
        a_obj.style.display   = "none";
        div_obj.style.display = "";
        frm_obj.pop_name.value  = a_obj.innerHTML;
    }
}

function annotate_population(id) {
    //
    // Fetch the marker annotation
    //
    var form_obj = document.getElementById(id + "_frm"); 
    var url = 
        form_obj.url.value + "&" + 
        "pop_name=" + escape(form_obj.pop_name.value);

    // 
    // Prepare and send XMLHttpRequest Object.
    //
    http_req = false;

    try {
	http_req = new XMLHttpRequest();
    } catch(e) {
	http_req = false;
    }

    if (http_req) {
	http_req.onreadystatechange = process_pop_annotation;
	http_req.open("GET", url, true);
	http_req.send("");
    }

    toggle_population(id);
}

function process_pop_annotation() {
    //
    // Possible readyState values:
    // 0 = uninitialized
    // 1 = loading
    // 2 = loaded
    // 3 = interactive
    // 4 = complete
    //
    if (http_req.readyState == 4) {

        // Check that the status is "OK"
        if (http_req.status == 200) {

            var xml_doc = http_req.responseXML;
            var id_obj   = xml_doc.getElementsByTagName("pop_id");
            var name_obj = xml_doc.getElementsByTagName("text");
            var pop_id   = id_obj[0].childNodes[0].nodeValue;
            var txt;
            if (name_obj[0].childNodes.length > 0)
                txt = name_obj[0].childNodes[0].nodeValue;
            else
                txt = "Population " + pop_id;

            var a_obj = document.getElementById(pop_id + "_pop");
            a_obj.innerHTML = txt;

        } else {
            alert("There was a problem retrieving the XML data:\n" + http_req.statusText);
        }
    }
}
