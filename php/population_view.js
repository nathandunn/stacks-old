var colors = {
    'a' : "#0074D9", // Blue
    'b' : "#FF4136", // Red 
    'c' : "#008000", // Dark Green
    'd' : "#FF851B", // Orange 
    'e' : "#001f3f", // Navy
    'f' : "#85144b", // Maroon 
    'g' : "#F012BE", // Fuchsia 
    'h' : "#39CCCC", // Teal 
    'i' : "#3D9970", // Olive 
    'j' : "#01FF70", // Lime 
    'k' : "#FFDC00", // Yellow 
    'l' : "#B10DC9", // Purple 
    'm' : "#111111", // Black 
    'n' : "#7FDBFF", // Aqua 
    'o' : "#AAAAAA", // Gray 
    'p' : "#DDDDDD", // Silver 
    'q' : "#2ECC40", // Green 
    'r' : "#c00000",
    's' : "#ffc600",
    't' : "#29356c",
    'u' : "#860000",
    'v' : "#dc6200",
    'w' : "#4b398e",
    'x' : "#008f56",
    'y' : "#bf1e25",
    'z' : "#4cb8ff",
    'A' : "#26081C",
    'B' : "#C44900",
    'C' : "#036016",
    'D' : "#936426",
    'E' : "#E87EA1",
    'F' : "#692E38",
    'G' : "#19423F",
    'H' : "#625F64",
    'I' : "#432534",
    'J' : "#323845",
    'K' : "#758E4F",
    'L' : "#B287A3",
    'M' : "#C0F8D1",
    'N' : "#03440C",
    'O' : "#8EDCE6",
    'P' : "#D5DCF9",
    'Q' : "#A7B0CA",
    'R' : "#0CBABA",
    'S' : "#725E54",
    'T' : "#EE2677",
    'U' : "#84A98C",
    'V' : "#FCCA46",
    'W' : "#A1C181",
    'X' : "#443627",
    'Y' : "#FE7F2D",
    'Z' : "#2EBE62"
};

marker_types = {
    'ab/--' : ['aa', 'bb', 'ab', '--'],
    '--/ab' : ['aa', 'bb', 'ab', '--'],
    'ab/aa' : ['aa', 'bb', 'ab', '--'],
    'aa/ab' : ['aa', 'bb', 'ab', '--'],
    'ab/ab' : ['aa', 'ab', 'bb', '--'],
    'ab/ac' : ['aa', 'ab', 'ac', 'bc', '--'],
    'ab/cd' : ['aa', 'bb', 'cc', 'dd', 'ac', 'ad', 'bc', 'bd', '--'],
    'aa/bb' : ['aa', 'bb', 'ab', '--'],
    'ab/cc' : ['aa', 'bb', 'ab', 'ac', 'bc', 'cc', '--'],
    'cc/ab' : ['aa', 'bb', 'ab', 'ac', 'bc', 'cc', '--']
};

function ajax_locus_population_view(id, root_path, url) {
    var tr_obj = document.getElementById(id); 
    var caret  = document.getElementById(id + "_img");

    if (tr_obj.style.display == "none") {
	tr_obj.style.display = "";
	caret.src            = root_path + "/caret-d.png"; 

	//
	// Make the AJAX query for JSON encoded data for this locus.
	//
	$.ajax({
	    dataType: "json",
	    url: url,
	    success: build_population_view
	});

    }
    else {
	tr_obj.style.display = "none";
	caret.src            = root_path + "/caret-u.png"; 
    }
}

function ajax_locus_sumstats_view(id, url) 
{
    var div_obj = document.getElementById(id + "_sumstat_div");

    if (div_obj.style.display == "none") {
	//
	// Make the AJAX query for JSON encoded data for this locus.
	//
	$.ajax({
	    dataType: "json",
	    url: url,
	    success: build_sumstats_view
	});
    }
}

function ajax_locus_fst_view(id, url) 
{
    var div_obj = document.getElementById(id + "_fst_div");

    if (div_obj.style.display == "none") {
	//
	// Make the AJAX query for JSON encoded data for this locus.
	//
	$.ajax({
	    dataType: "json",
	    url: url,
	    success: build_fst_view
	});
    }
}

function ajax_locus_hapstats_view(id, url) 
{
    var div_obj = document.getElementById(id + "_hapstat_div");

    if (div_obj.style.display == "none") {
	//
	// Make the AJAX query for JSON encoded data for this locus.
	//
	$.ajax({
	    dataType: "json",
	    url: url,
	    success: build_hapstats_view
	});
    }
}

function ajax_locus_phist_view(id, url) 
{
    var div_obj = document.getElementById(id + "_phist_div");

    if (div_obj.style.display == "none") {
	//
	// Make the AJAX query for JSON encoded data for this locus.
	//
	$.ajax({
	    dataType: "json",
	    url: url,
	    success: build_phist_view
	});
    }
}

function ajax_locus_stack_view(id, url) 
{
    var div_obj = document.getElementById(id + "_stacks_div");

    $("html, body").css("cursor", "wait");

    if (div_obj.style.display == "none") {
	//
	// Make the AJAX query for JSON encoded data for this locus.
	//
	$.ajax({
	    dataType: "json",
	    url: url,
	    success: build_stack_view
	});
    }
}

function build_population_view(json, status, jqXHR)
{
    var url = 
	json.path + "/stack_view.php" + 
	"?db="       + json.db + 
	"&batch_id=" + json.batch_id + 
	"&tag_id="   + json.id;

    var html = 
	"<div id=\"" + json.id + "_stacks_div\" class=\"stack_view\" style=\"display: none\"></div>\n" +
        "<div class=\"comp_view\">\n" +
	"<div id=\"" + json.id + "_viewstat_stacks\" class=\"viewstat top_viewstat\" style=\"display: none\" onclick=\"build_stack_list('" + json.id + "', '" + url + "')\">View Stacks</div>\n";

    var snps = write_snps_table(json);

    html += snps;

    //
    // Generate map of genotypes to haplotypes (for use in color picking).
    //
    var color_map = {};
    for(var i = 0; i < json.alleles.length; i++)
	color_map[json.alleles[i].htype] = json.alleles[i].gtype;

    var alleles = write_haplotypes_table(json);

    html +=
        alleles + 
	"</div>\n" +
	"<div class=\"matches\">\n";

    var tables = write_population_table(json, color_map);
    html += 
        tables +
        "</div>\n";

    $("#" + json.id + "_popview_div").html(html);
}

function write_snps_table(json) {
    var html = "";

    if (json.snps.length == 0)
	return html;

    html +=
	"<table id=\"snps_table_" + json.id + "\" class=\"snps\">\n" +
	"<tr>\n" +
	"<th colspan=\"3\">SNPs</th>\n" +
	"</tr>\n" +
	"<tr class=\"subhead\">\n" +
	"<th></th>\n" +
	"<th>Column</th>\n" +
	"<th>Alleles</th>\n" +
	"</tr>\n";

    for(var i = 0; i < json.snps.length; i++) {
	var snp = 
	    "<tr id=\"" + json.id + "_snps_tr_" + json.snps[i].col + 
	    "\" onmouseover=\"highlight_snp('" + json.id + "_" + json.snps[i].col + "')\" onmouseout=\"unhighlight_snp('" + json.id + "_" + json.snps[i].col + "')\">\n" +
	    "<td>" + (i + 1) + ".</td>" +
	    "<td>" + json.snps[i].col + "</td>" + 
	    "<td>" + json.snps[i].rank_1 + " / " +
	    json.snps[i].rank_2;
	if (json.snps[i].rank_3.length > 0)
	    snp += " / " + json.snps[i].rank_3;
	if (json.snps[i].rank_4.length > 0)
	    snp += " / " + json.snps[i].rank_4;
	snp += 
	    "<div id=\"" + json.id + "_snps_tr_" + json.snps[i].col + "_close\" class=\"close_x\" style=\"display: none\">x</div></td>" +
	    "</tr>";
	html += snp;
    }

    html += 
	"</table>\n";

    if (json.snpstat == 1) {
	url = json.path + "/sumstat_view.php?db=" + json.db + "&batch_id=" + json.batch_id + "&type=" + json.type + "&tag_id=" + json.id;
	html +=
	    "<div id=\"" + json.id + "_viewstat_sumstat\" class=\"viewstat\" onclick=\"ajax_locus_sumstats_view('" + json.id  + "', '" + url + "')\">Summary Stats</div>\n" +
	    "<div id=\"" + json.id + "_sumstat_div\" class=\"popup_stat\" style=\"display: none;\"></div>\n"; 
    }
    if (json.snpfst == 1) {
	url = json.path + "/fst_view.php?db=" + json.db + "&batch_id=" + json.batch_id + "&type=" + json.type + "&tag_id=" + json.id;
	html +=
	    "<div id=\"" + json.id + "_viewstat_fststat\" class=\"viewstat\" onclick=\"ajax_locus_fst_view('" + json.id  + "', '" + url + "')\">F<sub>ST</sub> Stats</div>\n" +
	    "<div id=\"" + json.id + "_fst_div\" class=\"popup_stat\" style=\"display: none;\"></div>\n";
    }

    return html;
}

function write_haplotypes_table(json) {
    var html = "";

    if (json.alleles.length == 0)
	return html;

    html += 
        "<table id=\"alleles_table_" + json.id + "\" class=\"alleles\">\n" +
	"<tr>\n" +
	"<th colspan=\"4\">Haplotypes</th>\n" +
	"</tr>\n" +
	"<tr class=\"subhead\">\n" +
	"<th></th>\n" +
	"<th>Genotype</th>\n" +
	"<th>Haplotype</th>\n" +
	"<th>Cnt</th>\n" +
	"</tr>\n";

    //
    // Count up the number of occurances of each haplotype.
    //
    haplotypes = {};
    for (var pop_id in json.populations) {
	for (var i = 0; i < json.populations[pop_id].length; i++)
	    for (var j = 0; j < json.populations[pop_id][i].obshap.length; j++)
		if (!(json.populations[pop_id][i].obshap[j].allele in haplotypes)) {
		    haplotypes[json.populations[pop_id][i].obshap[j].allele] = 1;
		} else {
		    haplotypes[json.populations[pop_id][i].obshap[j].allele]++;
		}
    }

    //
    // Write out the haplotype table.
    //
    for(var i = 0; i < json.alleles.length; i++) {
	html += 
	    "<tr id=\"" + json.id + "_alleles_tr_" + i + "\" onclick=\"highlight_haplotypes('" + json.id + "', '" + json.alleles[i].htype + "', '" + json.id + "_alleles_tr_" + i + "')\">\n" +
	    "<td>" + (i+1) + ".</td>\n" +
	    "<td style=\"color: " + colors[json.alleles[i].gtype] + ";\">" + json.alleles[i].gtype + "</td>" + 
	    "<td style=\"color: " + colors[json.alleles[i].gtype] + ";\"><div class=\"haplotype_def\">" + json.alleles[i].htype + "</div></td>" +
	    "<td>";
	if (json.alleles[i].htype in haplotypes)
	    html += haplotypes[json.alleles[i].htype];
	else
	    html += "0";
	html += 
	    "<div id=\"" + json.id + "_alleles_tr_" + i + "_close\" class=\"close_x\" style=\"display: none\">x</div></td>" +
	    "</tr>";
    }

    html +=
        "</table>\n";

    if (json.hapstat == 1) {
	url = json.path + "/hapstat_view.php?db=" + json.db + "&batch_id=" + json.batch_id + "&type=" + json.type + "&tag_id=" + json.id;
	html +=
	    "<div id=\"" + json.id + "_viewstat_hapstat\" class=\"viewstat\" onclick=\"ajax_locus_hapstats_view('" + json.id  + "', '" + url + "')\">Haplotype Stats</div>\n" +
	    "<div id=\"" + json.id + "_hapstat_div\" class=\"popup_stat\" style=\"display: none;\"></div>\n"; 
    }
    if (json.hapfst == 1) {
	url = json.path + "/phist_view.php?db=" + json.db + "&batch_id=" + json.batch_id + "&type=" + json.type + "&tag_id=" + json.id;
	html +=
	    "<div id=\"" + json.id + "_viewstat_phistat\" class=\"viewstat\" onclick=\"ajax_locus_phist_view('" + json.id  + "', '" + url + "')\">&Phi;<sub>ST</sub> Stats</div>\n" +
	    "<div id=\"" + json.id + "_phist_div\" class=\"popup_stat\" style=\"display: none;\"></div>\n";
    }

    return html;
}

function write_population_table(json, color_map) {
    var html = "";

    html +=
        "<input type=text id=\"" + json.id + "_selected\" value=\"0\" style=\"display: none;\" />" + 
        "<table id=\"locus_gtypes_" + json.id + "\" class=\"genotypes\">\n" +
	"<tr>\n" +
	"  <td colspan=\"10\" class=\"gtype_toggle\" style=\"text-align: right;\">\n" +
	"    <strong>View:</strong>\n" +
	"    <input type=\"checkbox\" id=\"hap_cb_" + json.id + "\" checked=\"checked\" onclick=\"toggle_genotypes(" + json.id + ", 'locus_gtypes_" + json.id + "', 'hap')\" />\n" +
	"    <a onclick=\"document.getElementById('hap_cb_" + json.id + "').click()\">Haplotypes</a>\n" +
	"    <input type=\"checkbox\" id=\"dep_cb_" + json.id + "\" onclick=\"toggle_genotypes(" + json.id + ", 'locus_gtypes_" + json.id + "', 'dep')\" /> \n" +
	"    <a onclick=\"document.getElementById('dep_cb_" + json.id + "').click()\">Allele Depths</a>\n" +
	"    <input type=\"checkbox\" id=\"lnl_cb_" + json.id + "\" onclick=\"toggle_genotypes(" + json.id + ", 'locus_gtypes_" + json.id + "', 'lnl')\" /> \n" +
	"    <a onclick=\"document.getElementById('lnl_cb_" + json.id + "').click()\">LnLs</a>\n";
    if (json.type == "map") {
	html += 
            "<input id=\"gen_cb_" + json.id + "\" type=\"checkbox\" onclick=\"toggle_genotypes(" + json.id + ", 'locus_gtypes_" + json.id + "', 'gen')\" /> " +
	    "<a onclick=\"document.getElementById('gen_cb_" + json.id + "').click()\">Genotypes</a>";
    }
    html +=
	"  </td>\n" +
	"</tr>\n";

    //
    // Determine the maximum number of columns in the table if less than 10.
    //
    var num_cols = 0;
    for (var pop_id in json.populations)
	num_cols = json.populations[pop_id].length > num_cols ? json.populations[pop_id].length : num_cols;
    num_cols = num_cols < 10 ? num_cols : 10;

    //
    // Iterate over each population.
    //
    for (var pop_id in json.populations) {
	html += "<tr>\n";

	if (Object.keys(json.populations).length > 1) {
	    html += 
	        "<td class=\"pop_id\" colspan=\"" + num_cols + "\">";
	    html += print_population_name(pop_id, json);
	    html += 
	        "</td>\n" +
		"</tr>\n" +
		"<tr>\n";
	}

	var col_index = 0;
	var hap_index = 0;

	for (var i = 0; i < json.populations[pop_id].length; i++) {
            col_index++;

	    var sample  = json.populations[pop_id][i];
	    var uniq_id = sample.sample_id + "|" + sample.obshap[0].tag_id;

	    html += 
	        "  <td id=\"" + uniq_id + "_td\">" +
		"<div class=\"title\" onclick=\"enqueue_sample('" + json.id + "', '" + uniq_id + "')\">" + 
		"<input type=\"checkbox\" name=\"stack_view\" id=\"" + uniq_id + "\" />" +
		sample['sample'].replace("_", " ") + "</div>\n";

	    var hap_strs = [];
	    var dep_strs = [];

	    for (var j = 0; j < sample.obshap.length; j++) {
		hap_index++;
		var c   = colors[color_map[sample.obshap[j].allele]];
		var url = 
		    json.path + "/stack_view.php" + 
		    "?db="       + json.db + 
		    "&batch_id=" + json.batch_id + 
		    "&samples="  + sample.sample_id + 
		    "&tags="     + sample.obshap[j].tag_id +
		    "&tag_id="   + json.id;
		var a =
		    "<a onclick=\"ajax_locus_stack_view('" + json.id + "', '" + url + "')\" " + 
		    "title=\"#" + sample.obshap[j].tag_id + "\" style=\"color: " + c + ";\">";
		hap_strs.push("<span id=\"haphi_" + sample.obshap[j].allele + "_" + hap_index + "\">" + a + sample.obshap[j].allele + "</a></span>");
		dep_strs.push("<span id=\"dephi_" + sample.obshap[j].allele + "_" + hap_index + "\">" + a + sample.obshap[j].depth  + "</a></span>");
	    }

	    var hap_str = hap_strs.join(" / ");
	    var dep_str = dep_strs.join(" / ");
	    var gen_str = "";

	    if (json.type == "map" && sample.genotype.length > 0) {

		var id      = "gtype_" + json.batch_id + "_" + json.id + "_" + sample.sample_id;
		var url     = json.path + "/correct_genotype.php?db=" + json.db + "&batch_id=" + json.batch_id + "&tag_id=" + json.id + "&sample_id=" + sample.sample_id;
		var jscript = "correct_genotype('" + id + "', '" + url + "')";
		var blur_js = "cancel_correction('" + id + "')";

  		var sel = generate_element_select(id, sample.marker, sample.genotype.toLowerCase(), jscript, blur_js);
		var genotype = sample.corrected == 1 ? 
   		    "<span class=\"corrected\">" + sample.genotype + "</span>" :
		    sample.genotype;

		gen_str = 
		    "<div id=\"gen_" + col_index + "\" style=\"display: none;\">" +
		    "<div id=\"" + id + "_div\"><a onclick=\"toggle_correction('" + id + "')\">" + genotype + "</a></div>\n" +
		    "  <div id=\""+ id + "_sel\" style=\"display: none;\">\n" +
		    sel +
		    "  </div>";
	    }

	    html +=
	        "<div class=\"haplotype\" id=\"hap_" + col_index + "\">" + hap_str + "</div>" +
		"<div id=\"dep_" + col_index + "\" style=\"display: none;\">" + dep_str + "</div>" + 
		"<div class=\"lnl\" id=\"lnl_" + col_index + "\" style=\"display: none;\">" + sample.lnl + "</div>" +
		gen_str +
		"</td>\n";

	    if (col_index % num_cols == 0)
		html += 
                    "</tr>\n" +
                    "<tr>\n";
	}
	
	while (col_index % num_cols != 0) {
	    html += "  <td></td>\n";
	    col_index++;
	}
    }

    html += "</tr>\n";

    return html;
}

function enqueue_sample(cat_id, id)
{
    var cb  = document.getElementById(id); 
    var td  = document.getElementById(id + "_td"); 
    var sel = document.getElementById(cat_id + "_selected"); 

    var checks = Number(sel.value);

    if (cb.checked == false) {
	cb.checked = true;
	td.style.backgroundColor = "#ffffa8";
	checks++;

    } else {
	cb.checked = false;
	td.style.backgroundColor = "#ffffff";
	checks--;
    }

    if (checks > 0)
	$("#" + cat_id + "_viewstat_stacks").css("display", "");
    else
	$("#" + cat_id + "_viewstat_stacks").css("display", "none");

    sel.value = checks;
}

function build_stack_list(id, url)
{
    var samples = new Array();
    var tags    = new Array();

    $('input[type=checkbox]').each(function () {
	if (this.name == "stack_view" && this.checked) {
	    var parts = this.id.split("|");
	    samples.push(parts[0]);
	    tags.push(parts[1]);
	}
   });

    url += "&samples=" + samples.join(",") + "&tags=" + tags.join(",");

    ajax_locus_stack_view(id, url);
}

function print_population_name(pop_id, json) {

    var pop_str = json.popkey[pop_id].length > 0 ? json.popkey[pop_id] : "Population " + pop_id;

    var html =
        "<div class=\"pop_annotation\" id=\"" + pop_id + "_pop\" onclick=\"toggle_population(" + pop_id + ")\">" + pop_str + "</div>\n" +
	"<div id=\"" + pop_id + "_div\" style=\"display: none; font-size: small;\">\n" +
	"  <form id=\"" + pop_id + "_frm\">\n" +
	"    <input type=\"hidden\" name=\"url\" value=\"" + json.path + "/annotate_population.php?db=" + json.db + "&batch_id=" + json.batch_id + "&pop_id=" + pop_id + "\" />\n" +
	"    <input type=\"input\" size=15 name=\"pop_name\" value=\"\" />\n" +
	"    <a onclick=\"annotate_population('" + pop_id + "')\">save</a>|<a onclick=\"toggle_population('" + pop_id + "')\">cancel</a>\n" +
	"  </form>\n" +
	"</div>\n";

    return html;
}

function generate_element_select(name, marker, selected_ele, change_js, blur_js) {
    var script_code = "";

    if (change_js.length > 0 && blur_js.length > 0) {
        script_code = " onchange=\"" + change_js + "\" onblur=\"" + blur_js + "\"";
    } else if (change_js.length > 0) {
        script_code = " onchange=\"" + change_js + "\"";
    } else if (blur_js.length > 0) {
        script_code = " onblur=\"" + blur_js + "\"";
    }

    var ctl = "  <select id=\"" + name + "\" name=\"" + name + "\"" + script_code + ">\n";

    for (var i = 0; i < marker_types[marker].length; i++) {
	var element = marker_types[marker][i];

        if (element == selected_ele) 
            ctl += "  <option selected=\"selected\">" + element + "</option>\n";
        else
            ctl += "  <option>" + element + "</option>\n";
    }

    ctl += "  </select>\n";

    return ctl;
}

function highlight_snp(id) {
    var span_obj = document.getElementById(id); 

    span_obj.className = "rank_1_hi";
}

function unhighlight_snp(id) {
    var span_obj = document.getElementById(id); 

    span_obj.className = "rank_1";
}

function unhighlight_haplotypes(cat_id) {
    $("#alleles_table_" + cat_id + " tr").css("border", "");
    $("#alleles_table_" + cat_id + " tr > td").css("background-color", "");
    $("#alleles_table_" + cat_id + " tr td:first-child").css("color", "");
    $("#alleles_table_" + cat_id + " tr td:last-child").css("color", "");

    var table_obj = document.getElementById("locus_gtypes_" + cat_id);

    var spans = table_obj.getElementsByTagName("span");
    for(var i = 0; i < spans.length; i++) {
            spans[i].className = "";
    }

    $("div.close_x").css("display", "none");
}

function highlight_haplotypes(cat_id, haplotype, tr_id) {
    unhighlight_haplotypes(cat_id)

    $("#" + tr_id).css("border", "4px solid #a93535");
    $("#" + tr_id + " > td").css("background-color", "#aaa");
    $("#" + tr_id + " td:first-child").css("color", "#fff");
    $("#" + tr_id + " td:last-child").css("color", "#fff");

    var table_obj = document.getElementById("locus_gtypes_" + cat_id);

    var spans = table_obj.getElementsByTagName("span");
    for(var i = 0; i < spans.length; i++) {
	var parts = spans[i].id.split("_");

        if (parts[0] == "haphi" &&  parts[1] == haplotype) { 
            spans[i].className = "haphi";
	}
    }

    $("#" + tr_id + "_close").css("display", "");
    $("#" + tr_id +  "_close").bind("click", {cat_id: cat_id}, function(event) {
	unhighlight_haplotypes(event.data.cat_id);
	event.stopPropagation();
    });
}

function build_sumstats_view(json, status, jqXHR)
{
    var html = 
	"<div id=\"" + json.id + "_sumstat_div_close\" class=\"popup_stat_close\" onclick=\"close_sumstat_view('" + json.id + "')\" title=\"Close window\">x</div>\n" +
        "<table id=\"sumstats_table_" + json.id + "\" class=\"alleles\">\n" +
	"<tr>\n" +
	"<th colspan=\"14\">SNP Summary Statistics</th>\n" +
	"</tr>\n" +
	"<tr class=\"subhead\">\n" +
	"  <th></th>\n" +
	"  <th>Pop</td>\n" +
	"  <th>BP</td>\n" +
	"  <th>Column</td>\n" +
	"  <th>Allele 1</td>\n" +
	"  <th>Allele 2</td>\n" +
	"  <th><acronym title=\"Major allele frequency\">P</acronym></td>\n" +
	"  <th><acronym title=\"Number of samples\">N</acronym></td>\n" +
	"  <th>Obs Het</td>\n" +
	"  <th>Obs Hom</td>\n" +
	"  <th>Exp Het</td>\n" +
	"  <th>Exp Hom</td>\n" +
	"  <th>&pi;</td>\n" +
	"  <th>F<sub>IS</sub></td>\n" +
	"</tr>\n";

    for (var i = 0; i < json.sumstats.length; i++) {
	html +=
	    "<tr onmouseover=\"highlight_snp('" + json.id + "_" + json.sumstats[i].col + "')\" onmouseout=\"unhighlight_snp('" + json.id + "_" + json.sumstats[i].col + "')\">\n" +
	    "<td>" + (i + 1) + ".</td>\n" +
	    "<td>" + json.sumstats[i].pop_id + "</td>\n" + 
	    "<td>" + json.sumstats[i].bp + "</td>\n" + 
	    "<td>" + json.sumstats[i].col + "</td>\n" + 
	    "<td>" + json.sumstats[i].p_allele + "</td>\n" + 
	    "<td>" + json.sumstats[i].q_allele + "</td>\n" + 
	    "<td>" + json.sumstats[i].p + "</td>\n" + 
	    "<td>" + json.sumstats[i].n + "</td>\n" + 
	    "<td>" + json.sumstats[i].obshet + "</td>\n" + 
	    "<td>" + json.sumstats[i].obshom + "</td>\n" + 
	    "<td>" + json.sumstats[i].exphet + "</td>\n" + 
	    "<td>" + json.sumstats[i].exphom + "</td>\n" + 
	    "<td>" + json.sumstats[i].pi + "</td>\n" + 
	    "<td>" + json.sumstats[i].fis + "</td>\n" + 
	    "</tr>\n";
    }

    html +=
        "</table>\n";

    $("#" + json.id + "_sumstat_div").html(html);

    //
    // We need to position the sumstat overlay. Identify the right-most 
    // border of the containing div.
    //
    parent_div = "#" + json.id + "_popview_div";
    t = $(parent_div + " table.snps").position();
    p = $(parent_div + " div.comp_view").position();
    top_coord   = t.top + 32;
    right_coord = p.left + $(parent_div + " div.comp_view").width() + 25;
    $("#" + json.id + "_sumstat_div").css({top: top_coord, left: right_coord});

    //
    // Set a maximum height for the containing div.
    //
    h = $(parent_div).height() - 50;
    $("#" + json.id + "_sumstat_div").css("max-height", h);

    //
    // Make sure the Fst and Hapstat divs are closed.
    //
    unhighlight_snp_row(json.id);
    close_hapstat_view(json.id);
    close_phist_view(json.id);
    close_stack_view(json.id);

    //
    // Bind the escape key to close this popup.
    //
    $(document).keyup(function(event){
	if(event.keyCode === 27)
            close_sumstat_view(json.id);
    });

    //
    // Display the Sumstats div.
    //
    $("#" + json.id + "_sumstat_div").css("display", "");
    $("#" + json.id + "_viewstat_sumstat").css({"color": "#a93535", "border-color": "#a93535", "background-color": "#ffffff"});
}

function close_sumstat_view(id)
{
    $("#" + id + "_sumstat_div").css("display", "none");
    $("#" + id + "_viewstat_sumstat").css({"color": "", "border-color": "", "background-color": ""});
}

function build_hapstats_view(json, status, jqXHR)
{
    var html = 
	"<div id=\"" + json.id + "_hapstat_div_close\" class=\"popup_stat_close\" onclick=\"close_hapstat_view('" + json.id + "')\" title=\"Close window\">x</div>\n" +
        "<table id=\"hapstats_table_" + json.id + "\" class=\"alleles\">\n" +
	"<tr>\n" +
	"<th colspan=\"6\">Haplotype Summary Statistics</th>\n" +
	"</tr>\n" +
	"<tr class=\"subhead\">\n" +
	"  <th></th>\n" +
	"  <th>Pop</td>\n" +
	"  <th>BP</td>\n" +
	"  <th><acronym title=\"Number of haplotypes sampled\">N</acronym></td>\n" +
	"  <th><acronym title=\"Number of distinct haplotypes\">Haplotype Cnt</acronym></td>\n" +
	"  <th>Gene Diversity</td>\n" +
	"  <th>Haplotype Diversity</td>\n" +
	"</tr>\n";

    for (var i = 0; i < json.hapstats.length; i++) {
	html +=
	    "<tr>\n" +
	    "<td>" + (i + 1) + ".</td>\n" +
	    "<td>" + json.hapstats[i].pop_id + "</td>\n" + 
	    "<td>" + json.hapstats[i].bp + "</td>\n" + 
	    "<td>" + json.hapstats[i].n + "</td>\n" + 
	    "<td>" + json.hapstats[i].hapcnt + "</td>\n" + 
	    "<td>" + json.hapstats[i].genediv + "</td>\n" + 
	    "<td>" + json.hapstats[i].hapdiv + "</td>\n" + 
	    "</tr>\n";
    }

    html +=
        "</table>\n";

    $("#" + json.id + "_hapstat_div").html(html);

    //
    // We need to position the sumstat overlay. Identify the right-most 
    // border of the containing div.
    //
    parent_div = "#" + json.id + "_popview_div";
    t = $(parent_div + " table.snps").position();
    p = $(parent_div + " div.comp_view").position();
    top_coord   = t.top + 32;
    right_coord = p.left + $(parent_div + " div.comp_view").width() + 25;
    $("#" + json.id + "_hapstat_div").css({top: top_coord, left: right_coord});

    //
    // Set a maximum height for the containing div.
    //
    h = $(parent_div).height() - 50;
    $("#" + json.id + "_hapstat_div").css("max-height", h);

    //
    // Make sure the Fst, Phist, and Sumstat divs are closed.
    //
    unhighlight_snp_row(json.id);
    close_sumstat_view(json.id);
    close_phist_view(json.id);
    close_stack_view(json.id);

    //
    // Bind the escape key to close this popup.
    //
    $(document).keyup(function(event){
	if(event.keyCode === 27)
            close_hapstat_view(json.id);
    });

    //
    // Display the Hapstats div.
    //
    $("#" + json.id + "_hapstat_div").css("display", "");
    $("#" + json.id + "_viewstat_hapstat").css({"color": "#a93535", "border-color": "#a93535", "background-color": "#ffffff"});
}

function close_hapstat_view(id)
{
    $("#" + id + "_hapstat_div").css("display", "none");
    $("#" + id + "_viewstat_hapstat").css({"color": "", "border-color": "", "background-color": ""});
}

function build_fst_view(json, status, jqXHR)
{
    var html = 
	"<div id=\"" + json.id + "_fst_div_close\" class=\"popup_stat_close\" onclick=\"unhighlight_snp_row('" + json.id + "')\" title=\"Close window\">x</div>\n" +
	"<table class=\"fst_key\">\n" +
	"<tr>\n" +
	"  <td class=\"fst\">F<sub>ST</sub></td>\n" +
	"  <td class=\"pi\">&pi;<sub>overall</sub></td>\n" +
	"</tr>\n" +
	"</table>\n";

    //
    // Create a population id => name map.
    //
    var pop_names = [];
    for (var pop_id in json.popkey)
	pop_names[pop_id] = json.popkey[pop_id];

    var keys    = Object.keys(pop_names);
    var pop_cnt = keys.length;
    var cols    = Object.keys(json.columns);

    for (var col in json.columns) {
	var matrix = [];

	for (var i = 0; i < keys.length; i++)
	    matrix[keys[i]] = [];

	for (var i = 0; i < json.columns[col].length; i++) {
	    var s = json.columns[col][i];
	    matrix[s.pid_1][s.pid_2] = s;
	    matrix[s.pid_2][s.pid_1] = s;
	}

	html +=
	    "<table id=\"" + json.id + "_" + col + "_fst\" class=\"fst\" style=\"display: none;\">\n" +
	    "<tr>\n" +
	    "  <td class=\"key\" style=\"background-color: #fff; border-left: 0px; border-right: 0px;\">&nbsp;</td>\n";

	for (var i = 0; i < pop_cnt; i++)
	    html += "  <td class=\"key\" style=\"min-width: " + (100 / (pop_cnt + 1) - 1) + "%\">" + pop_names[keys[i]] + "</td>\n";

	html += "</tr>\n";

	for (var i = 0; i < pop_cnt; i++) {
	    html +=
	        "<tr>\n" +
		"  <td class=\"pop_id\">" + pop_names[keys[i]] + "</td>\n";

	    var str = "";
	    var cssclass = "";

	    for (var j = 0; j < pop_cnt; j++) {
		if (i < j) {
		    cssclass = "class=\"pi\"";
		    if (keys[j] in matrix[keys[i]])
			str = matrix[keys[i]][keys[j]]['pi_o'];
		} else {
		    cssclass = "class=\"fst\"";
		    if (keys[j] in matrix[keys[i]])
			str = 
			    "<span title=\"" +
			    "Fisher's P value: " + matrix[keys[i]][keys[j]]['p'] + "; LOD: " + matrix[keys[i]][keys[j]]['lod'] + "\">" +
			    matrix[keys[i]][keys[j]]['fst'] + "</span>";
		}

		if (i == j)
		    html += "  <td class=\"diagonal\">&nbsp;</td>\n";
		else if (!(keys[j] in matrix[keys[i]]))
		    html += "  <td " + cssclass + ">&nbsp;</td>\n";
		else {
		    var id  = json.id + "_" + col + "_" + i + "_" + j;
		    var arg = "'" + json.id + "', '" + col + "', '" + i + "', '" + j + "'";
		    html +=
		        "  <td id=\"" + id + "\" " + cssclass + " onmouseover=\"highlight_fst_cells(" + arg + ")\" onmouseout=\"unhighlight_fst_cells(" + arg + ")\">" +
			str +
			"</td>\n";
		}
	    }

	    html += "</tr>\n";
	}

	html += "</table>\n";
    }

    html +=
        "</table>\n";

    $("#" + json.id + "_fst_div").html(html);

    //
    // We need to make the first Fst table visible.
    //
    highlight_snp_row(json.id, cols[0]);

    //
    // Hook up event handlers to handle clicking on other SNPs to change the Fst table in view.
    //
    for (var i = 0; i < cols.length; i++) {
	var tr_id = "#" + json.id + "_snps_tr_" + cols[i];
	$(tr_id).bind("click", {cat_id: json.id, col: cols[i]}, function(event) {
	    change_highlighted_snp_row(event.data.cat_id, event.data.col);
	});
    }

    //
    // We need to position the fst overlay. Identify the right-most 
    // border of the containing div.
    //
    parent_div = "#" + json.id + "_popview_div";
    t = $(parent_div + " table.snps").position();
    p = $(parent_div + " div.comp_view").position();
    top_coord   = t.top + 32;
    right_coord = p.left + $(parent_div + " div.comp_view").width() + 25;
    $("#" + json.id + "_fst_div").css({top: top_coord, left: right_coord});

    //
    // Set a maximum height for the containing div.
    //
    h = $(parent_div).height() - 32;
    $("#" + json.id + "_fst_div").css("max-height", h);

    //
    // Make sure the Sumstats, Hapstats, and Phist divs are closed.
    //
    close_sumstat_view(json.id);
    close_hapstat_view(json.id);
    close_phist_view(json.id);
    close_stack_view(json.id);

    //
    // Bind the escape key to close this popup.
    //
    $(document).keyup(function(event){
	if(event.keyCode === 27)
	    unhighlight_snp_row(json.id);
    });

    //
    // Display the Fst div.
    //
    $("#" + json.id + "_fst_div").css("display", "");
    $("#" + json.id + "_viewstat_fststat").css({"color": "#a93535", "border-color": "#a93535", "background-color": "#ffffff"});
}

function close_fst_view(id)
{
    $("#" + id + "_fst_div").css("display", "none");
    $("#" + id + "_viewstat_fststat").css({"color": "", "border-color": "", "background-color": ""});

    //
    // Remove handlers for clicking on SNP rows.
    //
    $("#snps_table_" + id + " tr").unbind("click");
}

function build_phist_view(json, status, jqXHR)
{
    var html = 
	"<div id=\"" + json.id + "_fst_div_close\" class=\"popup_stat_close\" onclick=\"close_phist_view('" + json.id + "')\" title=\"Close window\">x</div>\n" +
	"<table class=\"fst_key\">\n" +
	"<tr>\n" +
	"  <td class=\"fst\">&Phi;<sub>ST</sub></td>\n" +
	"  <td class=\"pi\">F<sub>ST</sub>&rsquo;</td>\n" +
	"</tr>\n" +
	"</table>\n";

    //
    // Create a population id => name map.
    //
    var pop_names = [];
    for (var pop_id in json.popkey)
	pop_names[pop_id] = json.popkey[pop_id];

    var keys    = Object.keys(pop_names);
    var pop_cnt = keys.length;

    var matrix = [];

    for (var i = 0; i < keys.length; i++)
	matrix[keys[i]] = [];

    for (var i = 0; i < json.phist.length; i++) {
	var s = json.phist[i];
	matrix[s.pid_1][s.pid_2] = s;
	matrix[s.pid_2][s.pid_1] = s;
    }

    html +=
        "<table id=\"" + json.id + "_phist\" class=\"fst\">\n" +
	"<tr>\n" +
	"  <td class=\"key\" style=\"background-color: #fff; border-left: 0px; border-right: 0px;\">&nbsp;</td>\n";

    for (var i = 0; i < pop_cnt; i++)
	html += "  <td class=\"key\" style=\"min-width: " + (100 / (pop_cnt + 1) - 1) + "%\">" + pop_names[keys[i]] + "</td>\n";

    html += "</tr>\n";

    for (var i = 0; i < pop_cnt; i++) {
	html +=
	"<tr>\n" +
	    "  <td class=\"pop_id\">" + pop_names[keys[i]] + "</td>\n";

	var str = "";
	var cssclass = "";

	for (var j = 0; j < pop_cnt; j++) {
	    if (i < j) {
		cssclass = "class=\"pi\"";
		if (keys[j] in matrix[keys[i]])
		    str = matrix[keys[i]][keys[j]]['fstp'];
	    } else {
		cssclass = "class=\"fst\"";
		if (keys[j] in matrix[keys[i]])
		    str = matrix[keys[i]][keys[j]]['phist'];
	    }

	    if (i == j)
		html += "  <td class=\"diagonal\">&nbsp;</td>\n";
	    else if (!(keys[j] in matrix[keys[i]]))
		html += "  <td " + cssclass + ">&nbsp;</td>\n";
	    else {
		var id  = json.id + "_" + i + "_" + j;
		var arg = "'" + json.id + "', '', '" + i + "', '" + j + "'";
		html +=
		    "  <td id=\"" + id + "\" " + cssclass + " onmouseover=\"highlight_fst_cells(" + arg + ")\" onmouseout=\"unhighlight_fst_cells(" + arg + ")\">" +
		    str +
		    "</td>\n";
	    }
	}

	html += "</tr>\n";
    }

    html += "</table>\n";

    $("#" + json.id + "_phist_div").html(html);

    //
    // We need to position the fst overlay. Identify the right-most 
    // border of the containing div.
    //
    parent_div = "#" + json.id + "_popview_div";
    t = $(parent_div + " table.snps").position();
    p = $(parent_div + " div.comp_view").position();
    top_coord   = t.top + 32;
    right_coord = p.left + $(parent_div + " div.comp_view").width() + 25;
    $("#" + json.id + "_phist_div").css({top: top_coord, left: right_coord});

    //
    // Set a maximum height for the containing div.
    //
    h = $(parent_div).height() - 32;
    $("#" + json.id + "_phist_div").css("max-height", h);

    //
    // Make sure the Sumstats, Hapstats, and Fst divs are closed.
    //
    close_sumstat_view(json.id);
    close_hapstat_view(json.id);
    unhighlight_snp_row(json.id);
    close_stack_view(json.id);

    //
    // Bind the escape key to close this popup.
    //
    $(document).keyup(function(event){
	if(event.keyCode === 27)
            close_phist_view(json.id);
    });

    //
    // Display the Fst div.
    //
    $("#" + json.id + "_phist_div").css("display", "");
    $("#" + json.id + "_viewstat_phistat").css({"color": "#a93535", "border-color": "#a93535", "background-color": "#ffffff"});
}

function close_phist_view(id)
{
    $("#" + id + "_phist_div").css("display", "none");
    $("#" + id + "_viewstat_phistat").css({"color": "", "border-color": "", "background-color": ""});
}

function highlight_fst_cells(id, snp, row, col) 
{
    var cell_1, cell_2;

    if (snp.length == 0) {
	cell_1 = document.getElementById(id + "_" + row + "_" + col);
	cell_2 = document.getElementById(id + "_" + col + "_" + row);
    } else {
	cell_1 = document.getElementById(id + "_" + snp + "_" + row + "_" + col);
	cell_2 = document.getElementById(id + "_" + snp + "_" + col + "_" + row);
    }

    if (row < col) {
	cell_1.style.backgroundColor = "#f1592a";
	cell_1.style.color = "#ffffff";
	cell_2.style.backgroundColor = "#24aae2";
	cell_2.style.color = "#ffffff";
    } else {
	cell_1.style.backgroundColor = "#24aae2";
	cell_1.style.color = "#ffffff";
	cell_2.style.backgroundColor = "#f1592a";
	cell_2.style.color = "#ffffff";
    }
}

function unhighlight_fst_cells(id, snp, row, col) 
{
    var cell_1, cell_2;
    if (snp.length == 0) {
	cell_1 = document.getElementById(id + "_" + row + "_" + col);
	cell_2 = document.getElementById(id + "_" + col + "_" + row);
    } else {
	cell_1 = document.getElementById(id + "_" + snp + "_" + row + "_" + col);
	cell_2 = document.getElementById(id + "_" + snp + "_" + col + "_" + row);
    }

    if (row < col) {
	cell_1.style.backgroundColor = "#fdc4b8";
	cell_1.style.color = "#000000";
	cell_2.style.backgroundColor = "#aee2f3";
	cell_2.style.color = "#000000";
    } else {
	cell_1.style.backgroundColor = "#aee2f3";
	cell_1.style.color = "#000000";
	cell_2.style.backgroundColor = "#fdc4b8";
	cell_2.style.color = "#000000";
    }
}

function unhighlight_snp_row(cat_id) {

    $("#snps_table_" + cat_id + " tr").css("border", "");
    $("#snps_table_" + cat_id + " tr > td").css("background-color", "");
    $("#snps_table_" + cat_id + " tr td:first-child").css("color", "");
    $("#snps_table_" + cat_id + " tr td:last-child").css("color", "");

    //
    // Make the Fst tables invisible.
    //
    close_fst_view(cat_id);

    $("div.close_x").css("display", "none");
}

function highlight_snp_row(cat_id, col) {
    unhighlight_snp_row(cat_id)

    var tr_id = cat_id + "_snps_tr_" + col;

    $("#" + tr_id).css("border", "4px solid #a93535");
    $("#" + tr_id + " > td").css("background-color", "#aaa");
    $("#" + tr_id + " td:first-child").css("color", "#fff");
    $("#" + tr_id + " td:last-child").css("color", "#fff");

    //
    // Make the Fst table visible.
    //
    $("#" + cat_id + "_" + col + "_fst").css("display", "");

    $("#" + tr_id + "_close").css("display", "");
    $("#" + tr_id +  "_close").bind("click", {cat_id: cat_id}, function(event) {
	unhighlight_snp_row(event.data.cat_id);
	event.stopPropagation();
    });
}

function change_highlighted_snp_row(cat_id, col) {
    $("#snps_table_" + cat_id + " tr").css("border", "");
    $("#snps_table_" + cat_id + " tr > td").css("background-color", "");
    $("#snps_table_" + cat_id + " tr td:first-child").css("color", "");
    $("#snps_table_" + cat_id + " tr td:last-child").css("color", "");

    $("div.close_x").css("display", "none");

    var tr_id = cat_id + "_snps_tr_" + col;

    $("#" + tr_id).css("border", "4px solid #a93535");
    $("#" + tr_id + " > td").css("background-color", "#aaa");
    $("#" + tr_id + " td:first-child").css("color", "#fff");
    $("#" + tr_id + " td:last-child").css("color", "#fff");

    //
    // Hide the Fst tables;
    //
    $("#" + cat_id + "_fst_div table").css("display", "none");

    //
    // Make the Fst table visible.
    //
    $("table.fst_key").css("display", "");
    $("#" + cat_id + "_" + col + "_fst").css("display", "");

    $("#" + tr_id + "_close").css("display", "");
    $("#" + tr_id +  "_close").bind("click", {cat_id: cat_id}, function(event) {
	unhighlight_snp_row(event.data.cat_id);
	event.stopPropagation();
    });
}

function build_stack_view(json, status, jqXHR)
{
    var html = 
	"<div id=\"" + json.id + "_stacks_div_close\" class=\"popup_stat_close\" onclick=\"close_stack_view('" + json.id + "')\" title=\"Close window\">x</div>\n";

    for (var i = 0; i < json.stacks.length; i++) {
	html +=
	    "<table class=\"stack\">\n" +
	    "<tr>\n" +
	    "  <th colspan=\"4\">" + json.stacks[i].sample_name + " <span style=\"font-size: smaller\">[#" + json.stacks[i].sample_id + "]</span>; Stack " + json.stacks[i].tag_id + "</th>\n" +
	    "</tr>\n" +
            "<tr>\n" +
	    "  <td class=\"num\">&nbsp;</td>\n" +
	    "  <td class=\"con\">&nbsp;</td>\n" +
	    "  <td class=\"id\">&nbsp;</td>\n" +
	    "  <td class=\"tag\">" + json.stacks[i].scale + "</td>\n" +
	    "</tr>\n" +
	    "<tr>\n" +
	    "  <td class=\"num\">&nbsp;</td>\n" +
	    "  <td class=\"con\">consensus</td>\n" +
	    "  <td class=\"id\"></td>\n" +
	    "  <td class=\"tag\">" + json.stacks[i].consensus + "</td>\n" +
	    "</tr>\n" +
	    "<tr>\n" +
	    "  <td class=\"num\">&nbsp;</td>\n" +
	    "  <td class=\"con\">model</td>\n" +
	    "  <td class=\"id\"></td>\n" +
	    "  <td class=\"tag\">" + json.stacks[i].model + "</td>\n" +
	    "</tr>\n";

	//
	// Print out the primary stack components.
	//
	seq_cnt = 1;
	for (var j = 0; j < json.stacks[i].primary.length; j++) {
	    var bg = (j + 1) % 2 == 1 ? "style=\"background-color: #dddddd;\"" : "";
	    for (var k = 0; k < json.stacks[i].primary[j].ids.length; k++) {
		html +=
	            "<tr>\n" +
		    "  <td class=\"num\">" + seq_cnt + ".</td>\n" +
		    "  <td class=\"primary\">primary</td>\n" +
		    "  <td class=\"id\">" + json.stacks[i].primary[j].ids[k] + "</td>\n" +
		    "  <td class=\"tag\" " + bg + ">" + json.stacks[i].primary[j].seq + "</td>\n" +
		    "</tr>\n";
		seq_cnt++;
	    }
	}

	//
	// Print out the secondary stack components.
	//
	for (var j = 0; j < json.stacks[i].secondary.length; j++) {
	    html +=
	        "<tr>\n" +
		"  <td class=\"num\">" + seq_cnt + ".</td>\n" +
		"  <td class=\"secondary\">secondary</td>\n" +
		"  <td class=\"id\">"  + json.stacks[i].secondary[j].id  + "</td>\n" +
		"  <td class=\"tag\">" + json.stacks[i].secondary[j].seq + "</td>\n" +
		"</tr>\n";
	    seq_cnt++;
	}

	html += "</table>\n";
    }

    $("#" + json.id + "_stacks_div").html(html);

    //
    // Set a maximum height/width for the containing div.
    //
    var parent_div = "#" + json.id + "_popview_div";
    var p = $(parent_div).position();
    var h = $(parent_div).parent().height();
    var w = $(window).height();
    h = h > w ? w : h;
    var div_h = Math.round(h * 0.96);
    $("#" + json.id + "_stacks_div").css("max-height", div_h);

    $("html, body").css("cursor", "auto");

    //
    // Display the Stacks div.
    //
    $("#" + json.id + "_stacks_div").css("display", "");

    //
    // Make sure the Fst, Phist, and Sumstat divs are closed.
    //
    unhighlight_snp_row(json.id);
    close_sumstat_view(json.id);
    close_hapstat_view(json.id);
    close_phist_view(json.id);

    //
    // Bind the escape key to close this popup.
    //
    $(document).keyup(function(event){
	if(event.keyCode === 27)
            close_stack_view(json.id);
    });
}

function close_stack_view(id)
{
    $("#" + id + "_stacks_div").css("display", "none");
}
