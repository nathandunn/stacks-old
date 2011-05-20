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

function toggle_export_popup(id) {
    var div_obj = document.getElementById(id);

    if (div_obj.style.display == "none") {
        div_obj.style.display = "";
    }
    else {
        div_obj.style.display = "none";
    }
}

function close_export_popup(id) {
    var div_obj = document.getElementById(id);

    div_obj.style.display = "none";
}

function export_data(id) {
    var form_obj = document.getElementById(id + "_frm"); 
    var otype, dtype;

    for (i = 0; i < form_obj.otype.length; i++)
        if (form_obj.otype[i].checked)
            otype = form_obj.otype[i].value;

    for (i = 0; i < form_obj.dtype.length; i++)
        if (form_obj.dtype[i].checked)
            dtype = form_obj.dtype[i].value;

    var url = 
        form_obj.url.value + "&" + 
        "email=" + form_obj.email.value + "&" +
        "dtype=" + dtype + "&" +
        "dlim="  + form_obj.dlim.value + "&" +
        "mtype=" + form_obj.mtype.value + "&" +
        "mcor="  + form_obj.mcor.value + "&" +
        "otype=" + otype;
    //alert(url);

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
	http_req.onreadystatechange = process_export;
	http_req.open("GET", url, true);
	http_req.send("");
    }

    var txt_obj = document.getElementById(id + "_txt"); 

    while (txt_obj.childNodes.length > 0)
        txt_obj.removeChild(txt_obj.firstChild); 
}

function gen_export_result(loci, email) {

    var div_obj = document.getElementById("export_popup_txt"); 
    var res_p   = document.createElement("p");

    res_p.innerHTML = "Exporting <span class=\"r\">" + loci + "</span> loci. " +
                      "E-mail will be sent to <span class=\"r\">" + email + "</span>" +
                      " when it is complete.";

    var app_obj = div_obj.appendChild(res_p);

    var res_p       = document.createElement("p");
    res_p.innerHTML = "<a onclick=\"close_export_popup('export_popup')\">close</a>";

    app_obj.appendChild(res_p);
}

function process_export() {
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

            var obj  = xml_doc.getElementsByTagName("loci");
            var loci = obj[0].childNodes[0].nodeValue;

            obj       = xml_doc.getElementsByTagName("email");
            var email = obj[0].childNodes[0].nodeValue;

            obj       = xml_doc.getElementsByTagName("msg");
            var msg   = obj[0].childNodes[0].nodeValue;

            //alert(msg);

            gen_export_result(loci, email);

        } else {
            alert("There was a problem retrieving the XML data:\n" + http_req.statusText);
        }
    }
}

function toggle_vis(form_id, name) {
    var form_obj = document.getElementById(form_id); 
    var g_obj    = document.getElementById('gopts'); 
    var h_obj    = document.getElementById('hopts'); 

    for(i = 0; i < form_obj.elements.length; i++)
        if (form_obj.elements[i].name  == name &&
	    form_obj.elements[i].value == "geno") 

            if (form_obj.elements[i].checked == true) {
		g_obj.style.display = "";
		h_obj.style.display = "none";
	    } else {
                g_obj.style.display = "none";
		h_obj.style.display = "";
	    }
}
